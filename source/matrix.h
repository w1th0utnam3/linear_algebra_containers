/*
	linear_algebra_containers/matrix header file
	Copyright (C) 2015  Fabian Löschner

	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License along
	with this program; if not, write to the Free Software Foundation, Inc.,
	51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef MATRIX
#define MATRIX

#include "matrixbase.h"

namespace lin_algebra {

/**
 * Matrix template with support of basic linear alegbra operations
 *
 * This template provides a basic matrix implementation for m x n matrices. It
 * is not intended for large scale computations with huge dimensions but instead
 * for smaller problems and calculations in computer graphics.
 * To access entries there is a matrix index operator implemented as the ()
 * function operator. This operator uses the standard mathematical row-major
 * way of adressing entries: mat(row,column). All operations are implemented
 * by using this row-major () operator. Therefore everything works like it is
 * defined in most linear algebra books.
 * Internally the elments of the matrix are stored in a column-major ordering in
 * a 1-dimensional data array. This storage is either directly accessible with
 * the implemented [] array subscript operator or via the data() method which
 * returns a pointer to the underliyng data.
 *
 * @tparam T Type used for the entries of the matrix. Must support basic
 * aritmethic operations.
 * @tparam m Number of Matrix::rows of the matrix.
 * @tparam n Number of columns of the matrix.
 */
template<typename T, size_t row_count_param, size_t column_count_param>
class Matrix : public MatrixBase<T,row_count_param,column_count_param>
{
private:
	typedef MatrixBase<T,Matrix::rows,Matrix::cols> MatrixBaseType;

public:
	//! The type of the matrix
	typedef Matrix<T,Matrix::rows,Matrix::cols> MatrixType;
	//! The type of the transposed matrix
	typedef Matrix<T,Matrix::cols,Matrix::rows> TransposedMatrixType;

	/**
	 * @brief Construct a matrix
	 *
	 * Constructs a matrix without initilization or a m*n parameter pack
	 * for initialization depending on the arguments of the constructor.
	 */
	template<typename ...Ts>
	Matrix(Ts... values)
		: MatrixBaseType(values...)
	{
	}

	/**
	 * @brief Create an identity matrix
	 *
	 * This method constructs a matrix with all entries set to zero except
	 * for the diagonal entries which are set to one.
	 * @return An identity matrix.
	 */
	static MatrixType createIdentity()
	{
		MatrixType result;
		result.toIdentity();
		return result;
	}

	/**
	 * @brief Set matrix to identity matrix
	 *
	 * Sets all entries to zero except for the diagonal entries which are
	 * set to one.
	 */
	MatrixType& toIdentity()
	{
		this->fill(T(0));

		const size_t smaller_dim = (Matrix::rows < Matrix::cols) ? Matrix::rows : Matrix::cols;
		for(size_t i = 0; i < smaller_dim; i++) {
			this->entries_[MatrixBaseType::index(i,i)] = 1;
		}

		return *this;
	}

	/**
	 * @brief Create transposed matrix
	 *
	 * Constructs the transposed of this matrix, i.e. creating matrix
	 * B of dimensions [n x m] from matrix A [m x n] with A*B = I where
	 * I is an identity matrix.
	 * @return The transposed of this matrix.
	 */
	TransposedMatrixType transposed() const
	{
		TransposedMatrixType result;
		for(size_t i = 0; i < Matrix::rows; i++) {
			for(size_t j = 0; j < Matrix::cols; j++) {
				result(j,i) = this->entries_[MatrixBaseType::index(i,j)];
			}
		}
		return result;
	}

	//! Adds the right matrix to the left.
	MatrixType operator+=(const MatrixType& rhs)
	{
		for(size_t i = 0; i < (Matrix::rows*Matrix::cols); i++) {
			this->entries_[i] += rhs.entries_[i];
		}
		return *this;
	}

	//! Substracts the right matrix from the left.
	MatrixType operator-=(const MatrixType& rhs)
	{
		for(size_t i = 0; i < (Matrix::rows*Matrix::cols); i++) {
			this->entries_[i] -= rhs.entries_[i];
		}
		return *this;
	}

	//! Scales the matrix by the specified factor.
	MatrixType operator*=(double factor)
	{
		for(size_t i = 0; i < (Matrix::rows*Matrix::cols); i++) {
			this->entries_[i] *= factor;
		}
		return *this;
	}

	//! Returns the matrix scaled by the specified factor.
	friend MatrixType operator*(const MatrixType& mat, double factor)
	{
		MatrixType result(mat);
		result *= factor;
		return result;
	}

	//! Returns the matrix scaled by the specified factor.
	friend MatrixType operator*(double factor, const MatrixType& mat)
	{
		MatrixType result(mat);
		result *= factor;
		return result;
	}

	//! Returns the sum of the two matrices. Matrix dimensions must agree.
	friend MatrixType operator+(const MatrixType& lhs, const MatrixType& rhs)
	{
		MatrixType result(lhs);
		result += rhs;
		return result;
	}

	//! Returns the difference of the two matrices. Matrix dimensions must agree.
	friend MatrixType operator-(const MatrixType& lhs, const MatrixType& rhs)
	{
		MatrixType result(lhs);
		result -= rhs;
		return result;
	}

	//! Returns the negated matrix.
	friend MatrixType operator-(const MatrixType& in)
	{
		return T(-1)*MatrixType(in);
	}
};

//! Returns the matrix product of two matrices. Matrix dimensions must agree. ([m x n]*[n x p] = [m x p])
template<typename T, size_t m, size_t n, size_t p>
inline Matrix<T,m,p> operator*(const Matrix<T,m,n>& lhs, const Matrix<T,n,p>& rhs)
{
	Matrix<T,m,p> result;
	result.zeros();

	for(size_t i = 0; i < m; i++) {
		for(size_t j = 0; j < p; j++) {
			for(size_t k = 0; k < n; k++) {
				result(i,j) += lhs(i,k) * rhs(k,j);
			}
		}
	}

	return result;
}

//! Simplified product when the matrix product is the same as a scalar product. ([1 x m]*[m x 1] = [1])
template<typename T, size_t m>
inline T operator*(const Matrix<T,1,m>& lhs, const Matrix<T,m,1>& rhs)
{
	T result(0);
	for(size_t i = 0; i < m; i++) {
		result += lhs(0,i) * rhs(i,0);
	}
	return result;
}

}

#endif // MATRIX
