///////////////////////////////////////////////////////////////////////////////
//
// CUDA function definitions
//
///////////////////////////////////////////////////////////////////////////////
#include "kernel.h"

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
/**
template <unsigned int blockSize>
__device__ void warpReduce(volatile int *sdata, unsigned int tid) {
if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
}
template <unsigned int blockSize>
__global__ void reduce6(int *g_idata, int *g_odata, unsigned int n) {
extern __shared__ int sdata[];
unsigned int tid = threadIdx.x;
unsigned int i = blockIdx.x*(blockSize * 2) + tid;
unsigned int gridSize = blockSize * 2 * gridDim.x;
sdata[tid] = 0;
while (i < n) { sdata[tid] += g_idata[i] + g_idata[i + blockSize]; i += gridSize; }
__syncthreads();
if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
if (blockSize >= 128) { if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }
if (tid < 32) warpReduce(sdata, tid);
if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}
*/
///////////////////////////////////////////////////////////////////////////////
// error check function
#define cuda_safe_call(ans) { simCudaAssert((ans), __FILE__, __LINE__); }
void simCudaAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess) {
		fprintf(stderr, "CUDA assert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

///////////////////////////////////////////////////////////////////////////////
// used namespace
using namespace simula;
using namespace simCuda;

///////////////////////////////////////////////////////////////////////////////
// local namespace
namespace {
	typedef MoleculeType::core_t kMoleculeType;
	kMoleculeType *tlist_d, *tlist_h;
	void copy_molecule_type()
	{
		simSize msize = molecules.type_num();
		tlist_h = (kMoleculeType*)malloc(msize * sizeof(kMoleculeType));
		for (simSize i = 0; i < msize; ++i) {
			tlist_h[i] = molecules.type(i).core();
			// deep copy dot pos
			simI3* dot_pos_d;
			simSize len = tlist_h[i].dot_num;
			cuda_safe_call(cudaMalloc((simI3**)&dot_pos_d, len * sizeof(simI3)));
			cuda_safe_call(cudaMemcpy(dot_pos_d, tlist_h[i].dot_pos, len * sizeof(simI3), cudaMemcpyHostToDevice));
			tlist_h[i].dot_pos = dot_pos_d;
		}
		cuda_safe_call(cudaMalloc((kMoleculeType**)&tlist_d, msize * sizeof(kMoleculeType)));
		cuda_safe_call(cudaMemcpy(tlist_d, tlist_h, msize * sizeof(kMoleculeType), cudaMemcpyHostToDevice));
	}
	/////////////////////////////////////////////////////////////////////////////
	// device variable
	kMolecule* mlist_d;
	kMolecule* mlist_h;

	/////////////////////////////////////////////////////////////////////////////
	// deep copy molecule list to device
	cudaError_t mlist_to_dev(simBool free_flag = false)
	{
		// allocate memory
		simSize msize = molecules.molecule_num();
		mlist_h = (kMolecule*)malloc(msize * sizeof(kMolecule));
		cuda_safe_call(cudaMalloc((void**)&mlist_d, msize * sizeof(kMolecule)));
		// deep copy molecule data into C struct
		for (simI1 i = 0; i < msize; ++i) {
			mlist_h[i].x = molecules.molecule(i).x();
			mlist_h[i].y = molecules.molecule(i).y();
			mlist_h[i].d = molecules.molecule(i).d();
			mlist_h[i].i = molecules.molecule(i).self_id();
			mlist_h[i].t = molecules.molecule(i).type_id();
		}
		// copy data into device
		cudaError_t err = cudaMemcpy(mlist_d, mlist_h, msize * sizeof(kMolecule), cudaMemcpyHostToDevice);
		// free mlist_h
		if (free_flag) {
			free(mlist_h);
		}
		return err;
	}

	/////////////////////////////////////////////////////////////////////////////
	// deep copy molecule list to host
	cudaError_t mlist_to_host()
	{
		// allocate memory
		simSize msize = molecules.molecule_num();
		// copy data back to host
		cudaError_t err = cudaMemcpy(mlist_h, mlist_d, msize * sizeof(kMolecule), cudaMemcpyDeviceToHost);
		// deep copy back to struct
		if (err = cudaSuccess) {
			for (simI1 i = 0; i < msize; ++i) {
				molecules.molecule(i).set_x(mlist_h[i].x);
				molecules.molecule(i).set_y(mlist_h[i].y);
				molecules.molecule(i).set_d(mlist_h[i].d);
			}
		}
		return err;
	}

	/////////////////////////////////////////////////////////////////////////////
	// kernel function
	__device__ __host__ simI1 pmod_k(simI1 x, simI1 n) {
		return ((x % n) + n) % n;
	}
	// ==> to check if its neighboring points are occupied
	__device__ simBool check_empty(kMoleculeType* tlist, simI1 sx, simI1 sy, kMolecule& ms, kMolecule& mt, simI1 xlen, simI1 ylen)
	{
		simBool empty = true;
		kMoleculeType& type_s = tlist[ms.t - 1];
		kMoleculeType& type_t = tlist[mt.t - 1];
		simSize len_s = type_s.dot_num;
		simSize len_t = type_t.dot_num;

		for (simI1 i = 0; i < len_s; ++i) {
			simI1 x_s, y_s;
			x_s = pmod_k(sx + type_s.dot_pos[i].x,xlen);
			y_s = pmod_k(sy + type_s.dot_pos[i].y,ylen);
			for (simI1 j = 0; j < len_t; ++j) {
				simI1 x_t, y_t;
				x_t = pmod_k(mt.x + type_t.dot_pos[j].x, xlen);
				y_t = pmod_k(mt.y + type_t.dot_pos[j].y, ylen);
				if (x_s == x_t && y_s == y_t) { empty = false; break; }
			}
		}

		return empty;
	}
	__global__ void addKernel(kMoleculeType* tlist, kMolecule* mlist, simI1* r, simI1 size, simI1 xsize, simI1 ysize)
	{
		simI1 idx = threadIdx.x;
		simI1 sx = mlist[idx].x, sy = mlist[idx].y;

		r[idx * 4 + 0] = idx+1;
		r[idx * 4 + 1] = idx+1;
		r[idx * 4 + 2] = idx+1;
		r[idx * 4 + 3] = idx+1;

		for (simI1 i = 0; i < size; ++i) {
			if (idx != i) {
				if (!check_empty(tlist, sx + 1, sy, mlist[idx], mlist[i], xsize, ysize)) { r[idx * 4 + 0] = 0; }
				if (!check_empty(tlist, sx, sy + 1, mlist[idx], mlist[i], xsize, ysize)) { r[idx * 4 + 1] = 0; }
				if (!check_empty(tlist, sx - 1, sy, mlist[idx], mlist[i], xsize, ysize)) { r[idx * 4 + 2] = 0; }
				if (!check_empty(tlist, sx, sy - 1, mlist[idx], mlist[i], xsize, ysize)) { r[idx * 4 + 3] = 0; }
			}
		}

	}

	simI1* result_d;
	simI1* result_h;

	// Helper function for using CUDA to add vectors in parallel.
	void funcCuda()
	{
		simI1 msize = molecules.molecule_num();

		// Choose which GPU to run on, change this on a multi-GPU system.
		cuda_safe_call(cudaSetDevice(0));

		// Allocate constant memory
		result_h = (simI1*)malloc(4 * msize * sizeof(simI1));
		cuda_safe_call(cudaMalloc((void**)&result_d, 4 * msize * sizeof(simI1)));

		// Copy data
		cuda_safe_call(mlist_to_dev(true));

		copy_molecule_type();

		// Launch a kernel on the GPU with one thread for each element.
		addKernel <<< 1, msize >>> (tlist_d, mlist_d, result_d, msize, sub.xlen(), sub.ylen());

		// Check for any errors launching the kernel
		cuda_safe_call(cudaGetLastError());

		// Check any errors encountered during the launch.
		cuda_safe_call(cudaDeviceSynchronize());

		// Copy output vector from GPU buffer to host memory.
		cuda_safe_call(cudaMemcpy(result_h, result_d, 4 * msize * sizeof(simI1), cudaMemcpyDeviceToHost));
	}

};

int simCuda::main_temp()
{
	// Test overlap.
	funcCuda();
	// print result
	for (simI1 i = 0; i < molecules.molecule_num(); ++i) {
		printf("{%d,%d,%d,%d}\n", result_h[4 * i + 0], result_h[4 * i + 1], result_h[4 * i + 2], result_h[4 * i + 3]);
	}

	// for tracing tools such as Nsight and Visual Profiler
	cuda_safe_call(cudaDeviceReset());
	return 0;
}