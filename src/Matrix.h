#pragma once
#include <complex>
#include <vector>
#include <numeric>
#include <cstring>
#include <assert.h>
#include <iostream>
#include "mkl.h"
#include "mkl_types.h"

// convert MKL_Complex16 to std::complex<double>
// convert MKL_Complex16 to std::complex<double>
inline std::complex<double> to_std_complex(const MKL_Complex16& c) {
    return std::complex<double>(c.real, c.imag);
}

inline MKL_Complex16 conjugate(const MKL_Complex16& c) {
    return MKL_Complex16({c.real, -c.imag});
}

inline bool operator==(const MKL_Complex16& base, const MKL_Complex16& other) {
    return (base.real == other.real && base.imag == other.imag);
}

inline MKL_Complex16& operator*=(MKL_Complex16& base, const MKL_Complex16& other) {
    base.real = base.real * other.real - base.imag * other.imag;
    base.imag = base.real * other.imag + base.imag * other.real;
    return base;
}

inline MKL_Complex16& operator+=(MKL_Complex16& base, const MKL_Complex16& other) {
    base.real += other.real;
    base.imag += other.imag;
    return base;
}

inline MKL_Complex16& operator-=(MKL_Complex16& base, const MKL_Complex16& other) {
    base.real -= other.real;
    base.imag -= other.imag;
    return base;
}

inline MKL_Complex16 operator+(const MKL_Complex16& base, const MKL_Complex16& other) {
    return {base.real + other.real, base.imag + other.imag};
}

inline MKL_Complex16 operator-(const MKL_Complex16& base, const MKL_Complex16& other) {
    return {base.real - other.real, base.imag - other.imag};
}

inline MKL_Complex16 operator-(const MKL_Complex16& base) {
    return {-base.real, -base.imag};
}

inline MKL_Complex16 operator*(const MKL_Complex16& base, const MKL_Complex16& other) {
    return {base.real * other.real - base.imag * other.imag,
                           base.real * other.imag + base.imag * other.real};
}

// template MKL_Complex operators

template<typename T>
MKL_Complex16 operator*(const MKL_Complex16& base, const T& other) {
    return {base.real * other, base.imag * other};
}

template<typename T>
class Matrix {

protected:

    T* data_ = nullptr;
    MKL_INT64 rows_;
    MKL_INT64 cols_;

public:

    Matrix() : data_(nullptr), rows_(0), cols_(0) {}
    Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
        if (data_ != nullptr) {
            mkl_free(data_);
        }
        data_ = (T*) mkl_malloc(rows * cols * sizeof(T), 64);
        for (int i = 0; i < rows * cols; ++i) {
            data_[i] = T(0);
        }
    }
    Matrix(const Matrix<T>& other) : rows_(other.rows_), cols_(other.cols_) {
        if (data_ != nullptr) {
            mkl_free(data_);
        }
        data_ = (T*) mkl_malloc(rows_ * cols_ * sizeof(T), 64);
		std::memcpy(data_, other.data_, rows_ * cols_ * sizeof(T));
    }

    ~Matrix() {
		if (data_ != nullptr) {
            mkl_free(data_);
        }
    }

    const T* data_ptr() {
        // for testing
        return data_;
    }

    void clear_matrix() {
        for (int i = 0; i < rows_ * cols_; ++i) {
            data_[i] = T(0);
        }
    }

    MKL_INT64 rows() const {
        return rows_;
    }

    MKL_INT64 cols() const {
        return cols_;
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
        // assignment operator
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

protected:

    // MKL_Complex16* data_;
    // MKL_INT64 rows_;
    // MKL_INT64 cols_;

    // data needed for inverse and determinant
    bool LU_decomposed_ = false; // flag to check if LU decomposition is done
    bool determinant_computed_ = false; // flag to check if determinant is computed
    T* LU_ = nullptr;
    MKL_INT* ipiv_ = nullptr;
    T determinant_ = {0, 0}; // determinant of the matrix

public:

    ComplexDoubleMatrix() : Matrix<T>() {}
    ComplexDoubleMatrix(int rows, int cols) {
        if (this->data_ != nullptr) {
            mkl_free(this->data_);
        }
        this->rows_ = rows;
        this->cols_ = cols;
        this->data_ = (T*) mkl_malloc(rows * cols * sizeof(T), 64);
        for (int i = 0; i < rows * cols; ++i) {
            this->data_[i] = T({0, 0});
        }
    }
    ComplexDoubleMatrix(const ComplexDoubleMatrix<T>& other) : Matrix<T>(other) {}

    ~ComplexDoubleMatrix() {
        if (LU_decomposed_) {
            mkl_free(LU_);
            mkl_free(ipiv_);
        }
    }

    void print_matrix(){
		// this is gross but I don't know how else to format print a complex number
		int WIDTH = 50;
		int token_width = 0;
		int num_spaces = 0;
		std::stringstream ss;
		for (int i = 0; i < this->rows_; ++i) {
			for (int j = 0; j < this->cols_-1; ++j) {
				ss.str("");
				ss << to_std_complex(this->data_[i * this->cols_ + j]);
				token_width = ss.str().length();
				num_spaces = WIDTH - token_width;
				if (num_spaces < 1) {
					num_spaces = 1;
				}
				std::cout << std::string(num_spaces, ' ') << to_std_complex(this->data_[i * this->cols_ + j]) << ",";
			}
			std::cout << to_std_complex(this->data_[i * this->cols_ + this->cols_ - 1]) << "\n";
		}
	}

    void clear_LU() {
        if (LU_decomposed_) {
            mkl_free(LU_);
            mkl_free(ipiv_);
            LU_decomposed_ = false;
            determinant_computed_ = false;
            determinant_ = {0, 0};
            LU_ = nullptr;
            ipiv_ = nullptr;
        }
    }

    void clear_matrix() {
        // clear the matrix
        for (int i = 0; i < this->rows_ * this->cols_; ++i) {
            this->data_[i] = T({0, 0});
        }
        clear_LU();
        
    }

    ComplexDoubleMatrix<T> get_slice(int row_start, int row_end, int col_start, int col_end) {
        // returns a slice of the matrix
        assert(row_start >= 0 && row_end <= this->rows_);
        assert(col_start >= 0 && col_end <= this->cols_);
        assert(row_start < row_end);
        assert(col_start < col_end);
        ComplexDoubleMatrix<T> result(row_end - row_start, col_end - col_start);
        for (int i = row_start; i < row_end; ++i) {
            for (int j = col_start; j < col_end; ++j) {
                result(i - row_start, j - col_start) = this->data_[i * this->cols_ + j];
            }
        }
        return result;
    }

    static ComplexDoubleMatrix<T> identity(int size) {
        // returns an identity matrix of given size
        ComplexDoubleMatrix<T> result(size, size);
        for (int i = 0; i < size; ++i) {
            result(i, i) = T({1.0, 0.0});
        }
        return result;
    }

    ComplexDoubleMatrix<T> get_conj_transpose() {
        // returns the conjugate transpose of the matrix
        ComplexDoubleMatrix<T> result(this->rows_, this->cols_);
        for (int i = 0; i < this->rows_; ++i) {
            for (int j = 0; j < this->cols_; ++j) {
                result(j, i) = conjugate(this->data_[i * this->cols_ + j]);
            }
        }
        return result;
    }

    // cblas routines

    void copy_column(const ComplexDoubleMatrix<T>& other, int src_col, int dest_col) {
        assert(src_col >= 0 && src_col < other.cols_);
        assert(dest_col >= 0 && dest_col < this->cols_);
        assert(other.rows_ == this->rows_);
        cblas_zcopy(other.rows_, &(other.data_[src_col]), other.cols_, &(this->data_[dest_col]), this->cols_);
    }

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

    void populate_LU() {
        // LU decomposition of the matrix. This is needed for the inverse and determinant
        // matrix should be square
        // we want to allocate LU_ and ipiv_ only once
        if (LU_decomposed_) {
            return;
        }
        LU_ = (T*) mkl_malloc(this->rows_ * this->cols_ * sizeof(T), 64);
        ipiv_ = (MKL_INT*) mkl_malloc(this->rows_ * this->cols_ * sizeof(MKL_INT), 64);
        std::memcpy(LU_, this->data_, this->rows_ * this->cols_ * sizeof(T));
        for (int i = 0; i < this->rows_ * this->cols_; ++i) {
            ipiv_[i] = 0;
        }
        MKL_INT info;
        info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, this->rows_, this->cols_, this->LU_, this->cols_, ipiv_);
        if (info != 0) {
            std::cerr << "Error in LU decomposition: " << info << std::endl;
            clear_LU();
            throw std::runtime_error("Error in LU decomposition");
        }
        this->LU_decomposed_ = true;
        this->compute_determinant();
    }

    T compute_inverse() {
        // computes the inverse of the matrix
        // matrix should be square
        // return the determinant since we'll get it for free
        if (!this->LU_decomposed_) {
            try {
                this->populate_LU();
            }
            catch( const std::runtime_error& e) {
                // if the matrix is singular, we cannot invert it
                std::cerr << "Matrix cannot be inverted: " << e.what() << std::endl;
                throw e;
            }
        }
        T determinant = this->determinant();
        MKL_INT info;
        std::memcpy(this->data_, this->LU_, this->rows_ * this->cols_ * sizeof(T));
        info = LAPACKE_zgetri(LAPACK_ROW_MAJOR, this->rows_, this->data_, this->cols_, this->ipiv_);
        clear_LU();
        if (info != 0) {
            std::cerr << "Matrix cannot be inverted: " << info << std::endl;
            throw std::runtime_error("Matrix cannot be inverted");
        }
        return determinant;
    }

    T compute_determinant() {
        // computes the determinant of the matrix
        // matrix should be square
        if (!LU_decomposed_) {
            populate_LU();
        }
        T result = {1.0, 0.0};
        for (int i = 0; i < this->rows_; ++i) {
            result *= this->LU_[i * this->cols_ + i];
        }
        this->determinant_ = result;
        this->determinant_computed_ = true;
        return result;
    }

    T determinant() {
        // returns the determinant of the matrix
        if (!determinant_computed_) {
            return compute_determinant();
        }
        return this->determinant_;
    }

    bool is_determinant_zero() {
        // checks if the determinant is zero
        // we need this because the absolute value of the determinant can be very small
        if (!determinant_computed_) {
            compute_determinant();
        }
        for (int i = 0; i < this->rows_; ++i) {
            if (std::abs(to_std_complex(this->data_[i*this->cols_ + i])) < 1e-10) {
                return true;
            }
        }
        return false;

    }

    T subdeterminant(std::vector<int>& rows, std::vector<int>& cols) {
        // computes the subdeterminant of the matrix
        // rows and cols should be of the same size
        assert(rows.size() == cols.size());
        assert(rows.size() <= this->rows_ && cols.size() <= this->cols_);
        ComplexDoubleMatrix<T> submatrix(rows.size(), cols.size());
        for (size_t i = 0; i < rows.size(); ++i) {
            for (size_t j = 0; j < cols.size(); ++j) {
                submatrix(i, j) = this->data_[rows[i] * this->cols_ + cols[j]];
            }
        }
        return submatrix.determinant();
    }

    // lapack routines

    std::pair<ComplexDoubleMatrix<T>, Matrix<double>> hermitian_diagonalize() {
        MKL_INT info;
        ComplexDoubleMatrix<T> result(this->rows_, this->cols_);
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

    // other operators

    ComplexDoubleMatrix<T> operator-(const ComplexDoubleMatrix<T>& other) {
        assert(this->rows_ == other.rows_);
        assert(this->cols_ == other.cols_);
        ComplexDoubleMatrix<T> result(this->rows_, this->cols_);
        for (int i = 0; i < this->rows_ * this->cols_; ++i) {
            result.data_[i] = this->data_[i] - other.data_[i];
        }
        return result;
    }
};