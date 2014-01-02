/********************************************************************************
 * File:   magUtils.c
 * Author: Brian J Smith <brian-j-smith@uiowa.edu>
 *
 * Created on June 21, 2010, 9:19 AM
 *
 * This file is part of the magma R package.
 *
 * magma is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * magma is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with magma.  If not, see <http://www.gnu.org/licenses/>.
 ********************************************************************************/

#include "magUtils.h"

#include <cublas.h>
#include <cuda.h>
#include <cuda_runtime_api.h>


void magCopyMatrix(int m, int n, double *a, int lda, double *b, int ldb)
{
   int i, j;

   for(i = 0; i < m; i++) {
      for(j = 0; j < n; j++) {
         a[i + j * lda] = b[i + j * ldb];
      }
   }
}


SEXP magGetGPU(SEXP a, SEXP b)
{
   SEXP gpu = GET_ATTR(a, install("gpu"));

   if(isNull(gpu)) gpu = GET_ATTR(b, install("gpu"));
   if(isNull(gpu)) error("must supply object of type \"magma\"");

   return gpu;
}


void magLoad()
{
   cublasInit();
   checkCublasError("initialization failed");
}


void magUnload()
{
   cublasShutdown();
}


int checkCudaError(const char *msg)
{
	cudaError_t err = cudaGetLastError();
	if(err != cudaSuccess)
		error("CUDA %s : %s\n", msg, cudaGetErrorString(err));
	return 0;
}

char *cublasGetErrorString(cublasStatus err)
{
	switch(err) {
		case CUBLAS_STATUS_SUCCESS :
			return "operation completed successfully";
		case CUBLAS_STATUS_NOT_INITIALIZED :
			return "CUBLAS library not initialized";
		case CUBLAS_STATUS_ALLOC_FAILED :
			return "resource allocation failed";
		case CUBLAS_STATUS_INVALID_VALUE :
			return "unsupported numerical value was passed to function";
		case CUBLAS_STATUS_ARCH_MISMATCH :
			return "function requires an architectural feature absent from \
			the architecture of the device";
		case CUBLAS_STATUS_MAPPING_ERROR :
			return "access to GPU memory space failed";
		case CUBLAS_STATUS_EXECUTION_FAILED :
			return "GPU program failed to execute";
		case CUBLAS_STATUS_INTERNAL_ERROR :
			return "an internal CUBLAS operation failed";
		default :
			return "unknown error type";
	}
}

int checkCublasError(const char *msg)
{
	cublasStatus err = cublasGetError();
	if(err != CUBLAS_STATUS_SUCCESS)
		error("CUBLAS %s : %s\n", msg, cublasGetErrorString(err));
	return 0;
}
