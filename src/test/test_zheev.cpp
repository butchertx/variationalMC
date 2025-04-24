#include "mkl.h"
#include "mkl_types.h"
#include <cstring>
#include <iostream>
#include <complex>
#include "Matrix.h"

int main(){
    MKL_INT info;
    MKL_INT64 N = 2;
    MKL_Complex16* data = (MKL_Complex16*) mkl_malloc(N * N * sizeof(MKL_Complex16), 64);
    MKL_Complex16* result = (MKL_Complex16*) mkl_malloc(N * N * sizeof(MKL_Complex16), 64);
    double * eigenvalues = (double*) mkl_malloc(N * sizeof(double), 64);

    data[0] = {0.0, 0.0};
    data[1] = {1.0, 0.0};
    data[2] = {1.0, 0.0};
    data[3] = {0.0, 0.0};

    std::memcpy(result, data, N * N * sizeof(MKL_Complex16));
    info = LAPACKE_zheev_64(LAPACK_ROW_MAJOR, 'V', 'U', N, result, N, eigenvalues);
    std::cout << "Eigenvalues: " << eigenvalues[0] << ", " << eigenvalues[1] << std::endl;

    ComplexDoubleMatrix<MKL_Complex16> data1(2, 2);
    data1(0, 0) = {0.0, 0.0};
    data1(0, 1) = {1.0, 0.0};
    data1(1, 0) = {1.0, 0.0};
    data1(1, 1) = {0.0, 0.0};
    auto eigensystem = data1.hermitian_diagonalize();
    std::cout << "Eigenvalues: " << eigensystem.second(0, 0) << ", " << eigensystem.second(1, 0) << std::endl;
    return 0;
}