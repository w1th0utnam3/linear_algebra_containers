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

#include <cstdlib>
#include <array>
#include <utility>

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
 * Internally the elments of the matrix are stored in a colum-major ordering in
 * a 1-dimensional data array. This storage is either directly accessible with
 * the implemented [] array subscript operator or via the data() method which
 * returns a pointer to the underliyng data.
 *
 * @tparam T Type used for the entries of the matrix. Must support basic
 * aritmethic operations.
 * @tparam m Number of rows of the matrix.
 * @tparam n Number of columns of the matrix.
 */
template<typename T, size_t row_count, size_t column_count>
class matrix
{

public:
	//! The number of rows of this matrix type.
	static const size_t m = row_count;
	//! The number of columns of this matrix type.
	static const size_t n = column_count;

	//! The type of the matrix
	typedef matrix<T,m,n> matrix_type;
	//! The type of the transposed matrix
	typedef matrix<T,n,m> transposed_matrix_type;

	/**
	 * @brief Construct an empty matrix
	 *
	 * All entries are uninitialized if T is of fundamental type or
	 * constructed using their default constructor if T is a complex
	 * type.
	 */
	matrix() = default;

	/**
	 * @brief Construct a matrix with specified entries
	 *
	 * This constructor initializes the matrix with the specified entries.
	 * It also allows using the braced initializer list syntax to initalize
	 * the matrix, e.g.: matrix<T,2,2> mat{{1,2,3,4}}.
	 */
	matrix(const std::array<T,n*m>& array)
		: entries(array)
	{
	}

	/**
	 * @brief Construct a matrix with specified entries
	 *
	 * This constructor initializes the matrix with the specified entries.
	 * It also allows using the braced initializer list syntax to initalize
	 * the matrix, e.g.: matrix<T,2,2> mat{{1,2,3,4}}.
	 * This constructor uses perfect forwarind to pass the an r-value to the
	 * constructor of the data array.
	 */
	matrix(std::array<T,n*m>&& array)
		: entries(std::forward<std::array<T,n*m>>(array))
	{
	}

	/**
	 * @brief Create an identity matrix
	 *
	 * This method constructs a matrix with all entries set to zero except
	 * for the diagonal entries which are set to one.
	 * @return Retuns an identity matrix.
	 */
	static matrix_type createIdentity()
	{
		matrix_type result;
		result.toIdentity();
		return result;
	}

	/**
	 * @brief 2d matrix index access operator
	 *
	 * This operator returns a reference to the entry of the matrix specified
	 * by the parameters. The operator uses row-major indexing.
	 * @param row The row of the entry.
	 * @param column The column of the entry.
	 * @return A reference to the specified entry.
	 */
	T& operator()(size_t row, size_t column)
	{
		return entries[row + column*m];
	}

	/**
	 * @brief 2d matrix index access operator
	 *
	 * This operator returns a const reference to the entry of the matrix
	 * specified by the arguments. The operator uses row-major indexing.
	 * @param row The row of the entry.
	 * @param column The column of the entry.
	 * @return A const reference to the specified entry.
	 */
	const T& operator()(size_t row, size_t column) const
	{
		return entries[row + column*m];
	}

	/**
	 * @brief 1d matrix index access operator
	 *
	 * This operator returns a reference to the i-th entry of the matrix. The
	 * operator uses column-major indexing i.e. in a 3x3 matrix the entries
	 * [0], [1] and [2] will be the first column of the matrix.
	 * @param index The one dimensional index of the entry.
	 * @return A reference to the specified entry.
	 */
	T& operator[](size_t i)
	{
		return entries[i];
	}

	/**
	 * @brief 1d matrix index access operator
	 *
	 * This operator returns a const reference to the i-th entry of the matrix.
	 * The operator uses column-major indexing i.e. in a 3x3 matrix the entries
	 * [0], [1] and [2] will be the first column of the matrix.
	 * @param index The one dimensional index of the entry.
	 * @return A const reference to the specified entry.
	 */
	const T& operator[](size_t i) const
	{
		return entries[i];
	}

	/**
	 * @brief Fill matrix with a value
	 *
	 * Sets value as the value for all the entries in the matrix.
	 * @param value The value that should be used to fill the matrix
	 */
	void fill(const T& value)
	{
		entries.fill(value);
	}

	/**
	 * @brief Fill matrix with zeros
	 *
	 * Sets all entries of the matrix to zero.
	 */
	void zeros()
	{
		entries.fill(T(0));
	}

	/**
	 * @brief Set matrix to identity matrix
	 *
	 * Sets all entries to zero except for the diagonal entries which are
	 * set to one.
	 */
	void toIdentity()
	{
		entries.fill(T(0));

		const size_t smaller_dim = (row_count < column_count) ? row_count : column_count;
		for(size_t i = 0; i < smaller_dim; i++) {
			(*this)(i,i) = 1;
		}
	}

	/**
	 * @brief Get pointer data
	 *
	 * Returns a pointer to the first element in the underlying data
	 * array which is used to store all entries of the matrix. The matrix
	 * uses a column-major memory layout therefore it is possble to pass
	 * this pointer to openGL methods.
	 * @return Pointer to the entries of the matrix.
	 */
	T* data()
	{
		return entries.data();
	}

	/**
	 * @brief Get pointer data
	 *
	 * Returns a const pointer to the first element in the underlying data
	 * array which is used to store all entries of the matrix. The matrix
	 * uses a column-major memory layout therefore it is possble to pass
	 * this pointer to openGL methods.
	 * @return Const pointer to the entries of the matrix.
	 */
	T* data() const
	{
		return entries.data();
	}

	/**
	 * @brief Create transposed matrix
	 *
	 * Constructs the transposed of this matrix, i.e. creating matrix
	 * B of dimensions [n x m] from matrix A [m x n] with A*B = I where
	 * I is an identity matrix.
	 * @return The transposed of this matrix.
	 */
	transposed_matrix_type transposed() const
	{
		transposed_matrix_type result;
		for(size_t i = 0; i < m; i++) {
			for(size_t j = 0; j < n; j++) {
				result(j,i) = (*this)(i,j);
			}
		}
		return result;
	}

	//! Adds the right matrix to the left.
	matrix_type operator+=(const matrix_type& rhs)
	{
		for(size_t i = 0; i < (m*n); i++) {
			entries[i] += rhs.entries[i];
		}
		return *this;
	}

	//! Substracts the right matrix from the left.
	matrix_type operator-=(const matrix_type& rhs)
	{
		for(size_t i = 0; i < (m*n); i++) {
			entries[i] -= rhs.entries[i];
		}
		return *this;
	}

	//! Scales the matrix by the specified factor.
	matrix_type operator*=(double factor)
	{
		for(size_t i = 0; i < (m*n); i++) {
			entries[i] *= factor;
		}
		return *this;
	}

	//! Returns whether two matrices have identical entries
	friend bool operator==(const matrix_type& lhs, const matrix_type& rhs)
	{
		for(size_t i = 0; i < n*m; i++) {
			if (lhs.entries[i] != rhs.entries[i]) return false;
		}
		return true;
	}

	//! Returns whether two matrices have not identical entries
	friend bool operator!=(const matrix_type& lhs, const matrix_type& rhs)
	{
		return !(lhs == rhs);
	}

	//! Returns the matrix scaled by the specified factor.
	friend matrix<T,m,n> operator*(const matrix<T,m,n>& mat, double factor)
	{
		matrix<T,m,n> result(mat);
		result *= factor;
		return result;
	}

	//! Returns the matrix scaled by the specified factor.
	friend matrix<T,m,n> operator*(double factor, const matrix<T,m,n>& mat)
	{
		matrix<T,m,n> result(mat);
		result *= factor;
		return result;
	}

	//! Returns the sum of the two matrices. Matrix dimensions must agree.
	friend matrix<T,m,n> operator+(const matrix<T,m,n>& lhs, const matrix<T,m,n>& rhs)
	{
		matrix<T,m,n> result(lhs);
		result += rhs;
		return result;
	}

	//! Returns the difference of the two matrices. Matrix dimensions must agree.
	friend matrix<T,m,n> operator-(const matrix<T,m,n>& lhs, const matrix<T,m,n>& rhs)
	{
		matrix<T,m,n> result(lhs);
		result -= rhs;
		return result;
	}

protected:
	std::array<T,m*n> entries;
};

//! Returns the matrix product of two matrices. Matrix dimensions must agree. ([m x n]*[n x p] = [m x p])
template<typename T, size_t m, size_t n, size_t p>
inline matrix<T,m,p> operator*(const matrix<T,m,n>& lhs, const matrix<T,n,p>& rhs)
{
	matrix<T,m,p> result;
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
inline T operator*(const matrix<T,1,m>& lhs, const matrix<T,m,1>& rhs)
{
	T result(0);
	for(size_t i = 0; i < m; i++) {
		result += lhs(0,i) * rhs(i,0);
	}
	return result;
}

#endif // MATRIX
