/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#include "build.h"
#include "assets/geometry/include/geometry.h"

#include "PreProcessor.h"

#undef USE_DOUBLE								// If enabled, double-precesion is used

#define DATA_DEGREE 0							// The order of the B-Spline used to splat in data for color interpolation
#define WEIGHT_DEGREE 2							// The order of the B-Spline used to splat in the weights for density estimation
#define NORMAL_DEGREE 2							// The order of the B-Spline used to splat in the normals for constructing the Laplacian constraints
#define DEFAULT_FEM_DEGREE 1					// The default finite-element degree
#define DEFAULT_FEM_BOUNDARY BOUNDARY_NEUMANN	// The default finite-element boundary type
#define DEFAULT_DIMENSION 3						// The dimension of the system

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "MyMiscellany.h"
#include "CmdLineParser.h"
#include "PPolynomial.h"
#include "FEMTree.h"
//#include "Ply.h"
#include "VertexFactory.h"
//#include "Image.h"
#include "RegularGrid.h"
#include "PoissonRecon.h"

struct PossionReconWrapper
{
	PoissonParams CommandLine;
	const PoissonInputStream* InputStream = nullptr;
	PoissonOutputStream* OutputStream = nullptr;

	//--

	const float DefaultPointWeightMultiplier = 2.f;

	double Weight( double v , double start , double end )
	{
		v = ( v - start ) / ( end - start );
		if     ( v<0 ) return 1.;
		else if( v>1 ) return 0.;
		else
		{
			// P(x) = a x^3 + b x^2 + c x + d
			//		P (0) = 1 , P (1) = 0 , P'(0) = 0 , P'(1) = 0
			// =>	d = 1 , a + b + c + d = 0 , c = 0 , 3a + 2b + c = 0
			// =>	c = 0 , d = 1 , a + b = -1 , 3a + 2b = 0
			// =>	a = 2 , b = -3 , c = 0 , d = 1
			// =>	P(x) = 2 x^3 - 3 x^2 + 1
			return 2. * v * v * v - 3. * v * v + 1.;
		}
	}

	template< unsigned int Dim, class Real >
	struct FEMTreeProfiler
	{
		double t;

		void start(void) { t = Time(), FEMTree< Dim, Real >::ResetLocalMemoryUsage(); }
		void print(const char* header) const
		{
		}
		void dumpOutput(const char* header) const
		{
		}
		void dumpOutput2(std::vector< std::string >& comments, const char* header) const
		{
		}
	};

	template< class Real , unsigned int Dim >
	XForm< Real , Dim+1 > GetBoundingBoxXForm( Point< Real , Dim > min , Point< Real , Dim > max , Real scaleFactor )
	{
		Point< Real , Dim > center = ( max + min ) / 2;
		Real scale = max[0] - min[0];
		for( int d=1 ; d<Dim ; d++ ) scale = std::max< Real >( scale , max[d]-min[d] );
		scale *= scaleFactor;
		for( int i=0 ; i<Dim ; i++ ) center[i] -= scale/2;
		XForm< Real , Dim+1 > tXForm = XForm< Real , Dim+1 >::Identity() , sXForm = XForm< Real , Dim+1 >::Identity();
		for( int i=0 ; i<Dim ; i++ ) sXForm(i,i) = (Real)(1./scale ) , tXForm(Dim,i) = -center[i];
		return sXForm * tXForm;
	}
	template< class Real , unsigned int Dim >
	XForm< Real , Dim+1 > GetBoundingBoxXForm( Point< Real , Dim > min , Point< Real , Dim > max , Real width , Real scaleFactor , int& depth )
	{
		// Get the target resolution (along the largest dimension)
		Real resolution = ( max[0]-min[0] ) / width;
		for( int d=1 ; d<Dim ; d++ ) resolution = std::max< Real >( resolution , ( max[d]-min[d] ) / width );
		resolution *= scaleFactor;
		depth = 0;
		while( (1<<depth)<resolution ) depth++;

		Point< Real , Dim > center = ( max + min ) / 2;
		Real scale = (1<<depth) * width;

		for( int i=0 ; i<Dim ; i++ ) center[i] -= scale/2;
		XForm< Real , Dim+1 > tXForm = XForm< Real , Dim+1 >::Identity() , sXForm = XForm< Real , Dim+1 >::Identity();
		for( int i=0 ; i<Dim ; i++ ) sXForm(i,i) = (Real)(1./scale ) , tXForm(Dim,i) = -center[i];
		return sXForm * tXForm;
	}

	template< typename Real , unsigned int Dim , typename AuxData >
	using InputOrientedPointStreamInfo = typename FEMTreeInitializer< Dim , Real >::template InputPointStream< VectorTypeUnion< Real , typename VertexFactory::NormalFactory< Real , Dim >::VertexType , AuxData > >;

	template< typename Real , unsigned int Dim , typename AuxData >
	using InputOrientedPointStream = typename InputOrientedPointStreamInfo< Real , Dim , AuxData >::StreamType;

	template< class Real , unsigned int Dim , typename AuxData >
	XForm< Real , Dim+1 > GetPointXForm( InputOrientedPointStream< Real , Dim , AuxData > &stream , Real width , Real scaleFactor , int& depth )
	{
		Point< Real , Dim > min , max;
		InputOrientedPointStreamInfo< Real , Dim , AuxData >::BoundingBox( stream , min , max );
		return GetBoundingBoxXForm( min , max , width , scaleFactor , depth );
	}
	template< class Real , unsigned int Dim , typename AuxData >
	XForm< Real , Dim+1 > GetPointXForm( InputOrientedPointStream< Real , Dim , AuxData > &stream , Real scaleFactor )
	{
		Point< Real , Dim > min , max;
		InputOrientedPointStreamInfo< Real , Dim , AuxData >::BoundingBox( stream , min , max );
		return GetBoundingBoxXForm( min , max , scaleFactor );
	}

	template< unsigned int Dim , typename Real >
	struct ConstraintDual
	{
		Real target , weight;
		ConstraintDual( Real t , Real w ) : target(t) , weight(w){ }
		CumulativeDerivativeValues< Real , Dim , 0 > operator()( const Point< Real , Dim >& p ) const { return CumulativeDerivativeValues< Real , Dim , 0 >( target*weight ); };
	};
	template< unsigned int Dim , typename Real >
	struct SystemDual
	{
		Real weight;
		SystemDual( Real w ) : weight(w){ }
		CumulativeDerivativeValues< Real , Dim , 0 > operator()( const Point< Real , Dim >& p , const CumulativeDerivativeValues< Real , Dim , 0 >& dValues ) const { return dValues * weight; };
		CumulativeDerivativeValues< double , Dim , 0 > operator()( const Point< Real , Dim >& p , const CumulativeDerivativeValues< double , Dim , 0 >& dValues ) const { return dValues * weight; };
	};
	template< unsigned int Dim >
	struct SystemDual< Dim , double >
	{
		typedef double Real;
		Real weight;
		SystemDual( Real w ) : weight(w){ }
		CumulativeDerivativeValues< Real , Dim , 0 > operator()( const Point< Real , Dim >& p , const CumulativeDerivativeValues< Real , Dim , 0 >& dValues ) const { return dValues * weight; };
	};

	template< typename VertexFactory, typename Index, class Real>
	void WritePolygons(const VertexFactory& vFactory, CoredMeshData< typename VertexFactory::VertexType, Index >* mesh, std::function< typename VertexFactory::VertexType(typename VertexFactory::VertexType) > xForm)
	{
		size_t nr_vertices = mesh->outOfCoreVertexNum() + mesh->inCoreVertices.size();
		size_t nr_faces = mesh->polygonNum();
		size_t nr_triangles = mesh->trianglePolygonNum();

		//OutputStream->positions.resize(nr_vertices * 3);
		//OutputStream->indices.reserve(nr_faces * 3);

		//OutputStream->geometry->allocateFaces(nr_triangles);
		//OutputStream->geometry->allocatePoints(nr_vertices);

		OutputStream->initializeMesh(nr_vertices, nr_triangles);

		//auto* vertexWritePtr = OutputStream->geometry->points().typedData();
		//auto* faceWritePtr = OutputStream->geometry->faces().typedData();

		mesh->resetIterator();

		std::vector< CoredVertexIndex< Index > > polygon;
		polygon.reserve(10);
		std::vector< int > indices;
		indices.reserve(10);

		for (size_t i = 0; i < mesh->inCoreVertices.size(); i++)
		{
			typename VertexFactory::VertexType vertex = xForm(mesh->inCoreVertices[i]);
			const float* pos = (const float*)&vertex;

			OutputStream->appendVertex(pos[1], pos[2], pos[3], 1.0f);// pos[0]);
			//OutputStream->appendVertex(coords[0], coords[1], coords[2], std::get< 0 >(vertex.data.psData).psData);

			//vertexWritePtr->v.x = pos[1];
			//vertexWritePtr->v.y = pos[2];
			//vertexWritePtr->v.z = pos[3];
			//++vertexWritePtr;
		}

		for (size_t i = 0; i < mesh->outOfCoreVertexNum(); i++)
		{
			typename VertexFactory::VertexType vertex = vFactory();
			mesh->nextOutOfCoreVertex(vertex);
			vertex = xForm(vertex);
			
			const float* pos = (const float*)&vertex;
			OutputStream->appendVertex(pos[1], pos[2], pos[3], 1.0f);// pos[0]);
			 
			//const auto coords = vertex.point.coords;
			//OutputStream->appendVertex(coords[0], coords[1], coords[2], std::get< 0 >(vertex.data.psData).psData);

			//vertexWritePtr->v.x = pos[1];
			//vertexWritePtr->v.y = pos[2];
			//vertexWritePtr->v.z = pos[3];
			//++vertexWritePtr;
		}

		for (size_t i = 0; i < nr_faces; i++)
		{
			polygon.clear();
			mesh->nextPolygon(polygon);

			indices.clear();
			indices.resize(polygon.size());

			for (int j = 0; j < polygon.size(); j++)
			{
				if (polygon[j].inCore) 
					indices[j] = polygon[j].idx;
				else                    
					indices[j] = polygon[j].idx + mesh->inCoreVertices.size();
			}

			for (int j = 2; j < indices.size(); ++j)
				OutputStream->appendTriangle(indices[0], indices[j - 1], indices[j]);
		}
	}

	template< typename Real , typename SetVertexFunction , typename InputSampleDataType , typename VertexFactory , unsigned int ... FEMSigs >
	void ExtractMesh
	(
		UIntPack< FEMSigs ... > ,
		FEMTree< sizeof ... ( FEMSigs ) , Real >& tree ,
		const DenseNodeData< Real , UIntPack< FEMSigs ... > >& solution ,
		Real isoValue ,
		const std::vector< typename FEMTree< sizeof ... ( FEMSigs ) , Real >::PointSample > *samples ,
		std::vector< InputSampleDataType > *sampleData ,
		const typename FEMTree< sizeof ... ( FEMSigs ) , Real >::template DensityEstimator< WEIGHT_DEGREE > *density ,
		const VertexFactory &vertexFactory ,
		const InputSampleDataType &zeroInputSampleDataType ,
		SetVertexFunction SetVertex ,
		std::vector< std::string > &comments ,
		XForm< Real , sizeof...(FEMSigs)+1 > unitCubeToModel
	)
	{
		static const int Dim = sizeof ... ( FEMSigs );
		typedef UIntPack< FEMSigs ... > Sigs;
		typedef typename VertexFactory::VertexType Vertex;

		static const unsigned int DataSig = FEMDegreeAndBType< DATA_DEGREE , BOUNDARY_FREE >::Signature;
		typedef typename FEMTree< Dim , Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;

		FEMTreeProfiler< Dim , Real > profiler;

	#if 0
		char tempHeader[1024];
		{
			char tempPath[1024];
			tempPath[0] = 0;
			if( TempDir.set ) strcpy( tempPath , TempDir.value );
			else SetTempDirectory( tempPath , sizeof(tempPath) );
			if( strlen(tempPath)==0 ) sprintf( tempPath , ".%c" , FileSeparator );
			if( tempPath[ strlen( tempPath )-1 ]==FileSeparator ) sprintf( tempHeader , "%sPR_" , tempPath );
			else                                                  sprintf( tempHeader , "%s%cPR_" , tempPath , FileSeparator );
		}
	#endif

		auto mesh = new CoredVectorMeshData< Vertex, node_index_type >();

		profiler.start();
		typename IsoSurfaceExtractor< Dim , Real , Vertex >::IsoStats isoStats;
		if( sampleData )
		{
			SparseNodeData< ProjectiveData< InputSampleDataType , Real > , IsotropicUIntPack< Dim , DataSig > > _sampleData = tree.template setExtrapolatedDataField< DataSig , false >( *samples , *sampleData , (DensityEstimator*)NULL );
			for( const RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type >* n = tree.tree().nextNode() ; n ; n=tree.tree().nextNode( n ) )
			{
				ProjectiveData< InputSampleDataType , Real >* clr = _sampleData( n );
				if( clr ) (*clr) *= (Real)pow( CommandLine.DataX.value , tree.depth( n ) );
			}
			isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< InputSampleDataType >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , density , &_sampleData , solution , isoValue , *mesh , zeroInputSampleDataType , SetVertex , !CommandLine.LinearFit.set , CommandLine.Normals.value==NORMALS_GRADIENTS , !CommandLine.NonManifold.set , CommandLine.PolygonMesh.set , false );
		}
	#if defined( __GNUC__ ) && __GNUC__ < 5
	#ifdef SHOW_WARNINGS
	#warning "you've got me gcc version<5"
	#endif // SHOW_WARNINGS
		else isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< InputSampleDataType >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , density , (SparseNodeData< ProjectiveData< InputSampleDataType , Real > , IsotropicUIntPack< Dim , DataSig > > *)NULL , solution , isoValue , *mesh , zeroInputSampleDataType , SetVertex , !LinearFit.set , CommandLine.Normals.value==NORMALS_GRADIENTS , !NonManifold.set , PolygonMesh.set , false );
	#else // !__GNUC__ || __GNUC__ >=5
		else isoStats = IsoSurfaceExtractor< Dim , Real , Vertex >::template Extract< InputSampleDataType >( Sigs() , UIntPack< WEIGHT_DEGREE >() , UIntPack< DataSig >() , tree , density , NULL , solution , isoValue , *mesh , zeroInputSampleDataType , SetVertex , !CommandLine.LinearFit.set , CommandLine.Normals.value==NORMALS_GRADIENTS , !CommandLine.NonManifold.set , CommandLine.PolygonMesh.set , false );
	#endif // __GNUC__ || __GNUC__ < 4
		//messageWriter( "Vertices / Polygons: %llu / %llu\n" , (unsigned long long)( mesh->outOfCoreVertexNum()+mesh->inCoreVertices.size() ) , (unsigned long long)mesh->polygonNum() );
		std::string isoStatsString = isoStats.toString() + std::string( "\n" );
		//messageWriter( isoStatsString.c_str() );
		if(CommandLine.PolygonMesh.set ) profiler.dumpOutput2( comments , "#         Got polygons:" );
		else                  profiler.dumpOutput2( comments , "#        Got triangles:" );

		std::vector< std::string > noComments;
		typename VertexFactory::Transform unitCubeToModelTransform( unitCubeToModel );

		WritePolygons< VertexFactory, node_index_type, Real >(vertexFactory, mesh, unitCubeToModelTransform);
		//PLY::WritePolygons< VertexFactory , node_index_type , Real , Dim >( Out.value , vertexFactory , mesh , ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE , NoComments.set ? noComments : comments , unitCubeToModelTransform );

		delete mesh;
	}

	template< typename InputSampleType >
	void ExtractData(std::vector< InputSampleType >& inCorePoints)
	{
		const char* posReadPtr = (const char*)InputStream->positions;
		const char* normalReadPtr = (const char*)InputStream->normals;

		inCorePoints.resize(this->InputStream->count);
		for (uint32_t i = 0; i < InputStream->count; ++i)
		{
			auto& pos = inCorePoints[i].template get<0>();
			pos[0] = ((float*)posReadPtr)[0];
			pos[1] = ((float*)posReadPtr)[1];
			pos[2] = ((float*)posReadPtr)[2];


			auto& normal = inCorePoints[i].template get<1>().template get<0>();
			normal[0] = ((float*)normalReadPtr)[0];
			normal[1] = ((float*)normalReadPtr)[1];
			normal[2] = ((float*)normalReadPtr)[2];

			posReadPtr += InputStream->positionStride;
			normalReadPtr += InputStream->normalStride;
		}
	}

	template< class Real , typename AuxDataFactory , unsigned int ... FEMSigs >
	void Execute( UIntPack< FEMSigs ... > , const AuxDataFactory &auxDataFactory )
	{
		static const int Dim = sizeof ... ( FEMSigs );
		typedef UIntPack< FEMSigs ... > Sigs;
		typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > Degrees;
		typedef UIntPack< FEMDegreeAndBType< NORMAL_DEGREE , DerivativeBoundary< FEMSignature< FEMSigs >::BType , 1 >::BType >::Signature ... > NormalSigs;
		static const unsigned int DataSig = FEMDegreeAndBType< DATA_DEGREE , BOUNDARY_FREE >::Signature;
		typedef typename FEMTree< Dim , Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;
		typedef typename FEMTree< Dim , Real >::template InterpolationInfo< Real , 0 > InterpolationInfo;
		using namespace VertexFactory;

		// The factory for constructing an input sample
		typedef Factory< Real , PositionFactory< Real , Dim > , Factory< Real , NormalFactory< Real , Dim > , AuxDataFactory > > InputSampleFactory;

		// The factory for constructing an input sample's data
		typedef Factory< Real , NormalFactory< Real , Dim > , AuxDataFactory > InputSampleDataFactory;

		// The input point stream information: First piece of data is the normal; the remainder is the auxiliary data
		typedef InputOrientedPointStreamInfo< Real , Dim , typename AuxDataFactory::VertexType > InputPointStreamInfo;

		// The type of the input sample
		typedef typename InputPointStreamInfo::PointAndDataType InputSampleType;

		// The type of the input sample's data
		typedef typename InputPointStreamInfo::DataType InputSampleDataType;

		typedef            InputDataStream< InputSampleType >  InputPointStream;
		typedef TransformedInputDataStream< InputSampleType > XInputPointStream;

		InputSampleFactory inputSampleFactory( PositionFactory< Real , Dim >() , InputSampleDataFactory( NormalFactory< Real , Dim >() , auxDataFactory ) );
		InputSampleDataFactory inputSampleDataFactory( NormalFactory< Real , Dim >() , auxDataFactory );

		typedef RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type > FEMTreeNode;
		typedef typename FEMTreeInitializer< Dim , Real >::GeometryNodeType GeometryNodeType;
		std::vector< std::string > comments;
		//messageWriter( comments , "*************************************************************\n" );
		//messageWriter( comments , "*************************************************************\n" );
		//messageWriter( comments , "** Running Screened Poisson Reconstruction (Version %s) **\n" , VERSION );
		//messageWriter( comments , "*************************************************************\n" );
		//messageWriter( comments , "*************************************************************\n" );
		//if( !CommandLine.Threads.set ) messageWriter( comments , "Running with %d threads\n" , CommandLine.Threads.value );

		bool needNormalData = CommandLine.DataX.value>0 && CommandLine.Normals.value;
		bool needAuxData = CommandLine.DataX.value>0 && auxDataFactory.bufferSize();

		XForm< Real, Dim + 1 > modelToUnitCube, unitCubeToModel;
		modelToUnitCube = XForm< Real , Dim+1 >::Identity();
		unitCubeToModel = XForm< Real, Dim + 1 >::Identity();

		double startTime = Time();
		Real isoValue = 0;

		FEMTree< Dim , Real > tree( MEMORY_ALLOCATOR_BLOCK_SIZE );
		FEMTreeProfiler< Dim , Real > profiler;

		size_t pointCount;

		Real pointWeightSum;
		std::vector< typename FEMTree< Dim , Real >::PointSample >* samples = new std::vector< typename FEMTree< Dim , Real >::PointSample >();
		DenseNodeData< GeometryNodeType , IsotropicUIntPack< Dim , FEMTrivialSignature > > geometryNodeDesignators;
		std::vector< InputSampleDataType >* sampleData = NULL;
		DensityEstimator* density = NULL;
		SparseNodeData< Point< Real , Dim > , NormalSigs >* normalInfo = NULL;
		Real targetValue = (Real)0.5;

		// Read in the samples (and color data)
		{
			profiler.start();

			sampleData = new std::vector< InputSampleDataType >();

			std::vector< InputSampleType > inCorePoints;
			ExtractData(inCorePoints);
			InputPointStream* pointStream = new MemoryInputDataStream< InputSampleType >(inCorePoints.size(), &inCorePoints[0]);

			typename InputSampleFactory::Transform _modelToUnitCube( modelToUnitCube );
			auto XFormFunctor = [&]( InputSampleType &p ){ p = _modelToUnitCube( p ); };
			XInputPointStream _pointStream( XFormFunctor , *pointStream );
			if( CommandLine.Width.value>0 )
			{
				modelToUnitCube = GetPointXForm< Real , Dim , typename AuxDataFactory::VertexType >( _pointStream , CommandLine.Width.value , (Real)( CommandLine.Scale.value>0 ? CommandLine.Scale.value : 1. ) , CommandLine.Depth.value ) * modelToUnitCube;
				if( !CommandLine.SolveDepth.set ) CommandLine.SolveDepth.value = CommandLine.Depth.value;
				if( CommandLine.SolveDepth.value>CommandLine.Depth.value )
				{
					WARN( "Solution depth cannot exceed system depth: " , CommandLine.SolveDepth.value , " <= " , CommandLine.Depth.value );
					CommandLine.SolveDepth.value = CommandLine.Depth.value;
				}
				if( CommandLine.FullDepth.value>CommandLine.Depth.value )
				{
					WARN( "Full depth cannot exceed system depth: " , CommandLine.FullDepth.value , " <= " , CommandLine.Depth.value );
					CommandLine.FullDepth.value = CommandLine.Depth.value;
				}
				if( CommandLine.BaseDepth.value>CommandLine.FullDepth.value )
				{
					if( CommandLine.BaseDepth.set ) WARN( "Base depth must be smaller than full depth: " , CommandLine.BaseDepth.value , " <= " , CommandLine.FullDepth.value );
					CommandLine.BaseDepth.value = CommandLine.FullDepth.value;
				}
			}
			else modelToUnitCube = CommandLine.Scale.value>0 ? GetPointXForm< Real , Dim , typename AuxDataFactory::VertexType >( _pointStream , (Real)CommandLine.Scale.value ) * modelToUnitCube : modelToUnitCube;

			{
				typename InputSampleFactory::Transform _modelToUnitCube( modelToUnitCube );
				auto XFormFunctor = [&]( InputSampleType &p ){ p = _modelToUnitCube( p ); };
				XInputPointStream _pointStream( XFormFunctor , *pointStream );
				auto ProcessDataWithConfidence = [&]( const Point< Real , Dim > &p , typename InputPointStreamInfo::DataType &d )
				{
					Real l = (Real)Length( d.template get<0>() );
					if( !l || !std::isfinite( l ) ) return (Real)-1.;
					return (Real)pow( l , CommandLine.Confidence.value );
				};
				auto ProcessData = []( const Point< Real , Dim > &p , typename InputPointStreamInfo::DataType &d )
				{
					Real l = (Real)Length( d.template get<0>() );
					if( !l || !std::isfinite( l ) ) return (Real)-1.;
					d.template get<0>() /= l;
					return (Real)1.;
				};
				if(CommandLine.Confidence.value>0 ) pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< InputSampleDataType >( tree.spaceRoot() , _pointStream , CommandLine.Depth.value , *samples , *sampleData , true , tree.nodeAllocators[0] , tree.initializer() , ProcessDataWithConfidence );
				else                     pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< InputSampleDataType >( tree.spaceRoot() , _pointStream , CommandLine.Depth.value , *samples , *sampleData , true , tree.nodeAllocators[0] , tree.initializer() , ProcessData );
			}

			unitCubeToModel = modelToUnitCube.inverse();
			delete pointStream;

			//messageWriter( "Input Points / Samples: %llu / %llu\n" , (unsigned long long)pointCount , (unsigned long long)samples->size() );
			profiler.dumpOutput2( comments , "# Read input into tree:" );
		}

		DenseNodeData< Real , Sigs > solution;
		{
			DenseNodeData< Real , Sigs > constraints;
			InterpolationInfo* iInfo = NULL;
			int solveDepth = CommandLine.Depth.value;

			tree.resetNodeIndices( 0 , std::make_tuple() );

			// Get the kernel density estimator
			{
				profiler.start();
				density = tree.template setDensityEstimator< 1 , WEIGHT_DEGREE >( *samples , CommandLine.KernelDepth.value , CommandLine.SamplesPerNode.value );
				profiler.dumpOutput2( comments , "#   Got kernel density:" );
			}

			// Transform the Hermite samples into a vector field
			{
				profiler.start();
				normalInfo = new SparseNodeData< Point< Real , Dim > , NormalSigs >();
				std::function< bool ( InputSampleDataType , Point< Real , Dim >& ) > ConversionFunction = []( InputSampleDataType in , Point< Real , Dim > &out )
				{
					Point< Real , Dim > n = in.template get<0>();
					Real l = (Real)Length( n );
					// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
					if( !l ) return false;
					out = n / l;
					return true;
				};
				std::function< bool ( InputSampleDataType , Point< Real , Dim >& , Real & ) > ConversionAndBiasFunction = [this]( InputSampleDataType in , Point< Real , Dim > &out , Real &bias )
				{
					Point< Real , Dim > n = in.template get<0>();
					Real l = (Real)Length( n );
					// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
					if( !l ) return false;
					out = n / l;
					bias = (Real)( log( l ) * CommandLine.ConfidenceBias.value / log( 1<<(Dim-1) ) );
					return true;
				};
	#if 1
				if(CommandLine.ConfidenceBias.value>0 ) *normalInfo = tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleData , density , CommandLine.BaseDepth.value , CommandLine.Depth.value , (Real)CommandLine.LowDepthCutOff.value , pointWeightSum , ConversionAndBiasFunction );
				else                         *normalInfo = tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleData , density , CommandLine.BaseDepth.value , CommandLine.Depth.value , (Real)CommandLine.LowDepthCutOff.value , pointWeightSum , ConversionFunction );
	#else
				if(CommandLine.ConfidenceBias.value>0 ) *normalInfo = tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleData , density , 0 , CommandLine.Depth.value , (Real)CommandLine.LowDepthCutOff.value , pointWeightSum , ConversionAndBiasFunction );
				else                         *normalInfo = tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleData , density , 0 , CommandLine.Depth.value , (Real)CommandLine.LowDepthCutOff.value , pointWeightSum , ConversionFunction );
	#endif
				ThreadPool::Parallel_for( 0 , normalInfo->size() , [&]( unsigned int , size_t i ){ (*normalInfo)[i] *= (Real)-1.; } );
				profiler.dumpOutput2( comments , "#     Got normal field:" );
				//messageWriter( "Point weight / Estimated Measure: %g / %g\n" , pointWeightSum , pointCount*pointWeightSum );
			}

			if( !CommandLine.Density.set ) delete density , density = NULL;
			if( !needNormalData && !needAuxData ) delete sampleData , sampleData = NULL;

			// Add the interpolation constraints
			if( CommandLine.PointWeight.value>0 )
			{
				profiler.start();
				if( CommandLine.ExactInterpolation.set ) iInfo = FEMTree< Dim , Real >::template       InitializeExactPointInterpolationInfo< Real , 0 > ( tree , *samples , ConstraintDual< Dim , Real >( targetValue , (Real)CommandLine.PointWeight.value * pointWeightSum ) , SystemDual< Dim , Real >( (Real)CommandLine.PointWeight.value * pointWeightSum ) , true , false );
				else                         iInfo = FEMTree< Dim , Real >::template InitializeApproximatePointInterpolationInfo< Real , 0 > ( tree , *samples , ConstraintDual< Dim , Real >( targetValue , (Real)CommandLine.PointWeight.value * pointWeightSum ) , SystemDual< Dim , Real >( (Real)CommandLine.PointWeight.value * pointWeightSum ) , true , 1 );
				profiler.dumpOutput2( comments , "#Initialized point interpolation constraints:" );
			}

			// Trim the tree and prepare for multigrid
			{
				profiler.start();
				constexpr int MAX_DEGREE = NORMAL_DEGREE > Degrees::Max() ? NORMAL_DEGREE : Degrees::Max();
				typename FEMTree< Dim , Real >::template HasNormalDataFunctor< NormalSigs > hasNormalDataFunctor( *normalInfo );
				auto hasDataFunctor = [&]( const FEMTreeNode *node ){ return hasNormalDataFunctor( node ); };
				if( geometryNodeDesignators.size() ) tree.template finalizeForMultigrid< MAX_DEGREE , Degrees::Max() >( CommandLine.BaseDepth.value , CommandLine.FullDepth.value , hasDataFunctor , [&]( const FEMTreeNode *node ){ return node->nodeData.nodeIndex<(node_index_type)geometryNodeDesignators.size() && geometryNodeDesignators[node]==GeometryNodeType::EXTERIOR; } , std::make_tuple( iInfo ) , std::make_tuple( normalInfo , density , &geometryNodeDesignators ) );
				else                                 tree.template finalizeForMultigrid< MAX_DEGREE , Degrees::Max() >( CommandLine.BaseDepth.value , CommandLine.FullDepth.value , hasDataFunctor , []( const FEMTreeNode * ){ return false; } , std::make_tuple( iInfo ) , std::make_tuple( normalInfo , density ) );

				profiler.dumpOutput2( comments , "#       Finalized tree:" );
			}
			// Add the FEM constraints
			{
				profiler.start();
				constraints = tree.initDenseNodeData( Sigs() );
				typename FEMIntegrator::template Constraint< Sigs , IsotropicUIntPack< Dim , 1 > , NormalSigs , IsotropicUIntPack< Dim , 0 > , Dim > F;
				unsigned int derivatives2[Dim];
				for( int d=0 ; d<Dim ; d++ ) derivatives2[d] = 0;
				typedef IsotropicUIntPack< Dim , 1 > Derivatives1;
				typedef IsotropicUIntPack< Dim , 0 > Derivatives2;
				for( int d=0 ; d<Dim ; d++ )
				{
					unsigned int derivatives1[Dim];
					for( int dd=0 ; dd<Dim ; dd++ ) derivatives1[dd] = dd==d ?  1 : 0;
					F.weights[d][ TensorDerivatives< Derivatives1 >::Index( derivatives1 ) ][ TensorDerivatives< Derivatives2 >::Index( derivatives2 ) ] = 1;
				}
				tree.addFEMConstraints( F , *normalInfo , constraints , solveDepth );
				profiler.dumpOutput2( comments , "#  Set FEM constraints:" );
			}

			// Free up the normal info
			delete normalInfo , normalInfo = NULL;

			// Add the interpolation constraints
			if( CommandLine.PointWeight.value>0 )
			{
				profiler.start();
				tree.addInterpolationConstraints( constraints , solveDepth , std::make_tuple( iInfo ) );
				profiler.dumpOutput2( comments , "#Set point constraints:" );
			}

			//messageWriter( "Leaf Nodes / Active Nodes / Ghost Nodes / Dirichlet Supported Nodes: %llu / %llu / %llu / %llu\n" , (unsigned long long)tree.leaves() , (unsigned long long)tree.nodes() , (unsigned long long)tree.ghostNodes() , (unsigned long long)tree.dirichletElements() );
			//messageWriter( "Memory Usage: %.3f MB\n" , float( MemoryInfo::Usage())/(1<<20) );
		
			// Solve the linear system
			{
				profiler.start();
				typename FEMTree< Dim , Real >::SolverInfo sInfo;
				sInfo.cgDepth = 0;
				sInfo.cascadic = true;
				sInfo.vCycles = 1;
				sInfo.iters = CommandLine.Iters.value;
				sInfo.cgAccuracy = CommandLine.CGSolverAccuracy.value;
				sInfo.verbose = CommandLine.Verbose.set;
				sInfo.showResidual = CommandLine.ShowResidual.set;
				sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE;
				sInfo.sliceBlockSize = 1;
				sInfo.baseVCycles = CommandLine.BaseVCycles.value;
				typename FEMIntegrator::template System< Sigs , IsotropicUIntPack< Dim , 1 > > F( { 0. , 1. } );
				solution = tree.solveSystem( Sigs() , F , constraints , CommandLine.SolveDepth.value , sInfo , std::make_tuple( iInfo ) );
				profiler.dumpOutput2( comments , "# Linear system solved:" );
				if( iInfo ) delete iInfo , iInfo = NULL;
			}
		}

		{
			profiler.start();
			double valueSum = 0 , weightSum = 0;
			typename FEMTree< Dim , Real >::template MultiThreadedEvaluator< Sigs , 0 > evaluator( &tree , solution );
			std::vector< double > valueSums( ThreadPool::NumThreads() , 0 ) , weightSums( ThreadPool::NumThreads() , 0 );
			ThreadPool::Parallel_for( 0 , samples->size() , [&]( unsigned int thread , size_t j )
			{
				ProjectiveData< Point< Real , Dim > , Real >& sample = (*samples)[j].sample;
				Real w = sample.weight;
				if( w>0 ) weightSums[thread] += w , valueSums[thread] += evaluator.values( sample.data / sample.weight , thread , (*samples)[j].node )[0] * w;
			}
			);
			for( size_t t=0 ; t<valueSums.size() ; t++ ) valueSum += valueSums[t] , weightSum += weightSums[t];
			isoValue = (Real)( valueSum / weightSum );
			if( !needNormalData && !needAuxData ) delete samples , samples = NULL;
			profiler.dumpOutput( "Got average:" );
			//messageWriter( "Iso-Value: %e = %g / %g\n" , isoValue , valueSum , weightSum );
		}
		
		//if( Out.set )
		{
			if( CommandLine.Normals.value )
			{
				if( CommandLine.Density.set )
				{
					typedef Factory< Real , PositionFactory< Real , Dim > , NormalFactory< Real , Dim > , ValueFactory< Real > , AuxDataFactory > VertexFactory;
					VertexFactory vertexFactory( PositionFactory< Real , Dim >() , NormalFactory< Real , Dim >() , ValueFactory< Real >() , auxDataFactory );
					if( CommandLine.Normals.value==1 )
					{
						auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p , v.template get<1>() = d.template get<0>() , v.template get<2>() = w , v.template get<3>() = d.template get<1>(); };
						ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
					}
					else if( CommandLine.Normals.value==2 )
					{
						auto SetVertex = [this]( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p , v.template get<1>() = -g/(1<<CommandLine.Depth.value) , v.template get<2>() = w , v.template get<3>() = d.template get<1>(); };
						ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
					}
				}
				else
				{
					typedef Factory< Real , PositionFactory< Real , Dim > , NormalFactory< Real , Dim > , AuxDataFactory > VertexFactory;
					VertexFactory vertexFactory( PositionFactory< Real , Dim >() , NormalFactory< Real , Dim >() , auxDataFactory );
					if( CommandLine.Normals.value==1 )
					{
						auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p                                                 , v.template get<1>() = d.template get<0>() , v.template get<2>() = d.template get<1>(); };
						ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
					}
					else if( CommandLine.Normals.value==2 )
					{
						auto SetVertex = [this]( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p                                                 , v.template get<1>() = -g/(1<<CommandLine.Depth.value) , v.template get<2>() = d.template get<1>(); };
						ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
					}
				}
			}
			else
			{
				if( CommandLine.Density.set )
				{
					typedef Factory< Real , PositionFactory< Real , Dim > , ValueFactory< Real > , AuxDataFactory > VertexFactory;
					VertexFactory vertexFactory( PositionFactory< Real , Dim >() , ValueFactory< Real >() , auxDataFactory );
					auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p , v.template get<1>() = w , v.template get<2>() = d.template get<1>(); };
					ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
				}
				else
				{
					typedef Factory< Real , PositionFactory< Real , Dim > , AuxDataFactory > VertexFactory;
					VertexFactory vertexFactory( PositionFactory< Real , Dim >() , auxDataFactory );
					auto SetVertex = []( typename VertexFactory::VertexType &v , Point< Real , Dim > p , Point< Real , Dim > g , Real w , InputSampleDataType d ){ v.template get<0>() = p , v.template get<1>() = d.template get<1>(); };
					ExtractMesh( UIntPack< FEMSigs ... >() , tree , solution , isoValue , samples , sampleData , density , vertexFactory , inputSampleDataFactory() , SetVertex , comments , unitCubeToModel );
				}
			}
			if( sampleData ){ delete sampleData ; sampleData = NULL; }
		}
		if( density ) delete density , density = NULL;
		//messageWriter( comments , "#          Total Solve: %9.1f (s), %9.1f (MB)\n" , Time()-startTime , FEMTree< Dim , Real >::MaxMemoryUsage() );
	}

	#ifndef FAST_COMPILE
	template< unsigned int Dim , class Real , BoundaryType BType , typename AuxDataFactory >
	void Execute( const AuxDataFactory &auxDataFactory )
	{
#if 0
		switch( Degree.value )
		{
			case 1: return Execute< Real >( IsotropicUIntPack< Dim , FEMDegreeAndBType< 1 , BType >::Signature >() , auxDataFactory );
			case 2: return Execute< Real >( IsotropicUIntPack< Dim , FEMDegreeAndBType< 2 , BType >::Signature >() , auxDataFactory );
	//		case 3: return Execute< Real >( IsotropicUIntPack< Dim , FEMDegreeAndBType< 3 , BType >::Signature >() , auxDataFactory );
	//		case 4: return Execute< Real >( IsotropicUIntPack< Dim , FEMDegreeAndBType< 4 , BType >::Signature >() , auxDataFactory );
			default: ERROR_OUT( "Only B-Splines of degree 1 - 2 are supported" );
		}
#endif
		return Execute< Real >(IsotropicUIntPack< Dim, FEMDegreeAndBType< 1, BType >::Signature >(), auxDataFactory);
	}

	template< unsigned int Dim , class Real , typename AuxDataFactory >
	void Execute( const AuxDataFactory &auxDataFactory )
	{
#if 0
		switch( BType.value )
		{
			case BOUNDARY_FREE+1:      return Execute< Dim , Real , BOUNDARY_FREE      >( auxDataFactory );
			case BOUNDARY_NEUMANN+1:   return Execute< Dim , Real , BOUNDARY_NEUMANN   >( auxDataFactory );
			case BOUNDARY_DIRICHLET+1: return Execute< Dim , Real , BOUNDARY_DIRICHLET >( auxDataFactory );
			default: ERROR_OUT( "Not a valid boundary type: " , BType.value );
		}
#endif
		Execute< Dim, Real, BOUNDARY_DIRICHLET      >(auxDataFactory);
	}
	#endif // !FAST_COMPILE

	int main()
	{
		if( CommandLine.MaxMemoryGB.value>0 ) 
			SetPeakMemoryMB( CommandLine.MaxMemoryGB.value<<10 );

		ThreadPool::DefaultChunkSize = CommandLine.ThreadChunkSize.value;
		ThreadPool::DefaultSchedule = (ThreadPool::ScheduleType)CommandLine.ScheduleType.value;
		ThreadPool::Init( (ThreadPool::ParallelType)CommandLine.ParallelType.value , CommandLine.Threads.value );

		if( !CommandLine.BaseDepth.set ) 
			CommandLine.BaseDepth.value = CommandLine.FullDepth.value;

		if( !CommandLine.SolveDepth.set ) 
			CommandLine.SolveDepth.value = CommandLine.Depth.value;

	#ifdef USE_DOUBLE
		typedef double Real;
	#else // !USE_DOUBLE
		typedef float  Real;
	#endif // USE_DOUBLE
	
		if( !CommandLine.PointWeight.set ) 
			CommandLine.PointWeight.value = DefaultPointWeightMultiplier * CommandLine.Degree.value;

		Execute< DEFAULT_DIMENSION, Real >(VertexFactory::EmptyFactory< Real >());

		ThreadPool::Terminate();
		return EXIT_SUCCESS;
	}
};

bool PoissonRecon(const PoissonParams& params, const PoissonInputStream& input, PoissonOutputStream* output)
{
	PossionReconWrapper wrapper;
	wrapper.CommandLine = params;
	wrapper.InputStream = &input;
	wrapper.OutputStream = output;

	try
	{
		wrapper.main();
	}
	catch (...)
	{
		return false;
	}

	return true;
}

PoissonParams::PoissonParams()
	: Degree(DEFAULT_FEM_DEGREE)
	, Depth(8)
	, KernelDepth(0)
	, SolveDepth(0)
	, Iters(8)
	, FullDepth(5)
	, BaseDepth(0)
	, BaseVCycles(1)
	, Normals(NORMALS_NONE)
	, BType(DEFAULT_FEM_BOUNDARY + 1)
	, MaxMemoryGB(0)
	, ParallelType((int)ThreadPool::THREAD_POOL)
	, ScheduleType((int)ThreadPool::DefaultSchedule)
	, ThreadChunkSize((int)ThreadPool::DefaultChunkSize)
	, Threads((int)std::thread::hardware_concurrency())
	, DataX(32.f)
	, SamplesPerNode(1.5f)
	, Scale(1.1f)
	, Width(0.f)
	, Confidence(0.f)
	, ConfidenceBias(0.f)
	, CGSolverAccuracy(1e-3f)
	, LowDepthCutOff(0.f)
{}