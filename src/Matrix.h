#pragma once
#include <complex>
#include <vector>
#include <numeric>
#include <cstring>
#include <assert.h>
#include <iostream>
#include "mkl.h"
#include "mkl_types.h"

bool operator==(const MKL_Complex16& base, const MKL_Complex16& other) {
    return (base.real == other.real && base.imag == other.imag);
}

MKL_Complex16& operator+=(MKL_Complex16& base, const MKL_Complex16& other) {
    base.real += other.real;
    base.imag += other.imag;
    return base;
}

MKL_Complex16& operator-=(MKL_Complex16& base, const MKL_Complex16& other) {
    base.real -= other.real;
    base.imag -= other.imag;
    return base;
}

template<typename T>
class Matrix {

protected:

    T* data_;
    MKL_INT64 rows_;
    MKL_INT64 cols_;

public:

    Matrix() : data_(nullptr), rows_(0), cols_(0) {}
    Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
        data_ = (T*) mkl_malloc(rows * cols * sizeof(T), 64);
        for (int i = 0; i < rows * cols; ++i) {
            data_[i] = T(0);
        }
    }
    Matrix(const Matrix<T>& other) : rows_(other.rows_), cols_(other.cols_) {
        data_ = (T*) mkl_malloc(rows_ * cols_ * sizeof(T), 64);
		std::memcpy(data_, other.data_, rows_ * cols_ * sizeof(T));
    }

    ~Matrix() {
		mkl_free(data_);
    }

    static void copy_row(Matrix<T>& dest, const Matrix<T>& src, int row_idx_dest, int row_idx_src, int column_boundary = -1) {
        assert(row_idx_src >= 0 && row_idx_src < src.rows_);
        assert(row_idx_dest >= 0 && row_idx_dest < dest.rows_);
        assert(src.cols_ >= column_boundary);
        if (column_boundary == -1) {
            column_boundary = src.cols_;
        }
        assert(dest.cols_ >= column_boundary);
        std::memcpy(&(dest.data_[row_idx_dest * dest.cols_]), &(src.data_[row_idx_src * src.cols_]), dest.cols_ * sizeof(T));
    }

    T& operator()(int row, int col) {
        assert(row >= 0 && row < rows_);
        assert(col >= 0 && col < cols_);
        return data_[row * cols_ + col];
    }

    bool operator==(const Matrix<T>& other) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            return false;
        }
        for (int i = 0; i < rows_ * cols_; ++i) {
            if (data_[i] == other.data_[i]) {
                continue;
            }
            else {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const Matrix<T>& other) const {
        return !(*this == other);
    }
};

template<typename T>
class ComplexDoubleMatrix : public Matrix<T> {

public:

    ComplexDoubleMatrix() : Matrix<T>() {}
    ComplexDoubleMatrix(int rows, int cols) {
        this->rows_ = rows;
        this->cols_ = cols;
        this->data_ = (T*) mkl_malloc(rows * cols * sizeof(T), 64);
        for (int i = 0; i < rows * cols; ++i) {
            this->data_[i] = T({0, 0});
        }
    }
    ComplexDoubleMatrix(const ComplexDoubleMatrix<T>& other) : Matrix<T>(other) {}

    ~ComplexDoubleMatrix() { /* Destructor will automatically call the base class destructor */ }

    // cblas routines

    Matrix<double> vector_norm(int row_start, int row_end) {
        // computes the vector norm of each column, but only for a range of rows
        // this is used to compute the occupation number of single particle orbitals
        // row_end is exclusive
        assert(row_start >= 0 && row_end <= this->rows_);
        assert(row_start < row_end);
        Matrix<double> result(1, this->cols_);
        for (int i = 0; i < this->cols_; ++i) {
            result(0, i) = cblas_dznrm2(row_end - row_start, &(this->data_[row_start * this->cols_ + i]), this->cols_);
        }
        return result;
    }

    Matrix<double> vector_norm() {
        // computes the vector norm of each column
        // returns a row vector of the norms
        return this->vector_norm(0, this->rows_);
    }

    // lapack routines

    std::pair<ComplexDoubleMatrix<T>, Matrix<double>> hermitian_diagonalize() {
        MKL_INT info;
        ComplexDoubleMatrix<T> result(*this);
        Matrix<double> eigenvalues(this->rows_, 1);
        std::memcpy(result.data_, this->data_, this->rows_ * this->cols_ * sizeof(T));
        info = LAPACKE_zheev_64(LAPACK_ROW_MAJOR, 'V', 'U', result.rows_, result.data_, result.rows_, eigenvalues.data_);
        if (info != 0) {
            std::cerr << "Error in diagonalization: " << info << std::endl;
        }
        return std::pair<ComplexDoubleMatrix<T>, Matrix<double>>(result, eigenvalues);
    }

    ComplexDoubleMatrix<T> operator*(const ComplexDoubleMatrix<T>& other) {
        assert(this->cols_ == other.rows_);
        ComplexDoubleMatrix<T> result(this->rows_, other.cols_);
        MKL_Complex16 alpha = { 1.0, 0.0 }, beta = { 0.0, 0.0 };
        cblas_zgemm3m_64(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->rows_, other.cols_, this->cols_,
                    &alpha, this->data_, this->cols_, other.data_, other.cols_, &beta, result.data_, result.cols_);
        return result;
    }
};