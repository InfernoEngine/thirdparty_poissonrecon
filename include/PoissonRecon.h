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

#pragma once
#include "CmdLineParser.h"

enum NormalType
{
	NORMALS_NONE,
	NORMALS_SAMPLES,
	NORMALS_GRADIENTS,
	NORMALS_COUNT
};

struct PoissonInputStream
{
	const void* positions = nullptr;
	uint32_t positionStride = 0;

	const void* normals = nullptr;
	uint32_t normalStride = 0;

	float pointOffset = 0.0f;

	uint32_t count = 0;
};

struct PoissonOutputStream
{
	virtual void initializeMesh(uint32_t numVertices, uint32_t numTriangles) = 0;
	virtual void appendVertex(float x, float y, float z, float d) = 0;
	virtual void appendTriangle(int a, int b, int c) = 0;
};

struct PoissonParams
{
	cmdLineReadable ShowResidual;
	cmdLineReadable NoComments;
	cmdLineReadable PolygonMesh;
	cmdLineReadable NonManifold;
	cmdLineReadable ASCII;
	cmdLineReadable Density;
	cmdLineReadable LinearFit;
	cmdLineReadable PrimalGrid;
	cmdLineReadable ExactInterpolation;
	cmdLineReadable Colors;
	cmdLineReadable InCore;
	cmdLineReadable NoDirichletErode;
	cmdLineReadable Verbose;

	cmdLineParameter< int > Degree;
	cmdLineParameter< int > Depth;
	cmdLineParameter< int > KernelDepth;
	cmdLineParameter< int > SolveDepth;
	cmdLineParameter< int > EnvelopeDepth;
	cmdLineParameter< int > Iters;
	cmdLineParameter< int > FullDepth;
	cmdLineParameter< int > BaseDepth;
	cmdLineParameter< int > BaseVCycles;
	cmdLineParameter< int > Normals;
	cmdLineParameter< int > BType;
	cmdLineParameter< int > MaxMemoryGB;
	cmdLineParameter< int > ParallelType;
	cmdLineParameter< int > ScheduleType;
	cmdLineParameter< int > ThreadChunkSize;
	cmdLineParameter< int > Threads;

	cmdLineParameter< float > DataX;
	cmdLineParameter< float > SamplesPerNode;
	cmdLineParameter< float > Scale;
	cmdLineParameter< float > Width;
	cmdLineParameter< float > Confidence;
	cmdLineParameter< float > ConfidenceBias;
	cmdLineParameter< float > CGSolverAccuracy;
	cmdLineParameter< float > LowDepthCutOff;
	cmdLineParameter< float > PointWeight;

	PoissonParams();
};

#if !defined(POISSON_RECON_SHARED)
	#define POISSON_RECON_API
#elif defined(POISSON_RECON_SHARED_EXPORTS)
	#define POISSON_RECON_API __declspec(dllexport)
#else
	#define POISSON_RECON_API __declspec(dllimport)
#endif

extern POISSON_RECON_API bool PoissonRecon(const PoissonParams& params, const PoissonInputStream& input, PoissonOutputStream* output);

