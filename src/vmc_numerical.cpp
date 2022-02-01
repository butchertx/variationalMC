#include "vmc_numerical.h"
#include <assert.h>
#include <algorithm>

std::vector<double> solve_linear_system(std::vector<std::vector<double>> A, std::vector<double> b) {
	//solve the linear system A*x = b.  Return x.
	//assume the major index of A is rows
	assert(A.size() == b.size());
	double* Ajk_l, * bk_l;
	lapack_int n, nrhs, lda, * ipiv, ldb;
	int M = A.size(),  N = A[0].size();

	Ajk_l = (double*)mkl_malloc(M * N * sizeof(double), 64);
	bk_l = (double*)mkl_malloc(M * sizeof(double), 64);
	ipiv = (int*)mkl_malloc(std::min({ M, N }) * sizeof(int), 64);

	for (int k = 0; k < N; ++k) {
		for (int j = 0; j < M; ++j) {
			Ajk_l[j * N + k] = A[j][k];
		}
		bk_l[k] = b[k];
	}

	//first get LU factorization
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, M, N, Ajk_l, N, ipiv);

	//then solve
	LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', M, 1, Ajk_l, N, ipiv, bk_l, 1);

	//transfer results to a std::vector
	std::vector<double> result(M);
	for (int i = 0; i < M; ++i) {
		result[i] = bk_l[i];
	}

	mkl_free(Ajk_l);
	mkl_free(bk_l);
	mkl_free(ipiv);

	return result;
}