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
#include <initializer_list>

/**
 * Template matrix type with support of basic matrix operations.
 *
 * To access entries there is a matrix index operator implemented as the ()
 * function operator. This operator uses the standard mathematical row-major
 * way of adressing entries: mat(row,column). All operations are implemented
 * by using this row-major () operator. Therefore everything works like it is
 * defined in most linear algebra scripts.
 * Internally the elments of the matrix are stored in a colum-major ordering in
 * a 1-dimensional data array. This storage is either accessible via the
 * implemented [] array subscript operator or via the data() method which
 * returns a pointer to the underliyng data.
 *
 * @tparam T Type used for the entries in the matrix.
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
	//! The type of the transposed
	typedef matrix<T,n,m> transposed_matrix_type;

	//! Constructs a matrix with non-initialized entries.
	matrix() = default;

	//! Constructs a matrix with the specified entries.
	matrix(std::array<T,n*m> array)
		: entries(array)
	{
	}

	//! Constructs a matrix initalized to an identity matrix
	inline static matrix_type createIdentity()
	{
		matrix_type result;
		result.eye();
		return result;
	}

	//! Returns a reference to the specified entry.
	inline T& operator()(size_t row, size_t column)
	{
		return entries[row + column*m];
	}

	//! Returns a const reference to the specified entry.
	inline const T& operator()(size_t row, size_t column) const
	{
		return entries[row + column*m];
	}

	//! Returns a reference to the specified entry.
	inline T& operator[](size_t i)
	{
		return entries[i];
	}

	//! Returns a const reference to the specified entry.
	inline const T& operator[](size_t i) const
	{
		return entries[i];
	}

	//! Fills the matrix with the specified value.
	inline void fill(T&& value)
	{
		entries.fill(std::forward<T>(value));
	}

	//! Fills the matrix with zeros.
	inline void zeros()
	{
		entries.fill(T(0));
	}

	//! Sets the matrix to an identity matrix.
	inline void eye()
	{
		entries.fill(T(0));

		const size_t smaller_dim = (row_count < column_count) ? row_count : column_count;
		for(size_t i = 0; i < smaller_dim; i++) {
			(*this)(i,i) = 1;
		}
	}

	//! Returns a pointer to the data of the matrix
	inline T* data()
	{
		return entries.data();
	}

	//! Returns a const pointer to the data of the matrix
	inline T* data() const
	{
		return entries.data();
	}

	//! Returns the transposed of the matrix.
	inline transposed_matrix_type transposed() const
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
	inline matrix_type operator+=(const matrix_type& rhs)
	{
		for(size_t i = 0; i < (m*n); i++) {
			entries[i] += rhs.entries[i];
		}
		return *this;
	}

	//! Substracts the right matrix from the left.
	inline matrix_type operator-=(const matrix_type& rhs)
	{
		for(size_t i = 0; i < (m*n); i++) {
			entries[i] -= rhs.entries[i];
		}
		return *this;
	}

	//! Scales the matrix by the specified factor.
	inline matrix_type operator*=(double factor)
	{
		for(size_t i = 0; i < (m*n); i++) {
			entries[i] *= factor;
		}
		return *this;
	}

	//! Returns whether two matrices have identical entries
	inline friend bool operator==(const matrix_type& lhs, const matrix_type& rhs)
	{
		for(size_t i = 0; i < n*m; i++) {
			if (lhs.entries[i] != rhs.entries[i]) return false;
		}
		return true;
	}

	//! Returns whether two matrices have not identical entries
	inline friend bool operator!=(const matrix_type& lhs, const matrix_type& rhs)
	{
		return !(lhs == rhs);
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

//! Simplified product when the matrix product is the same as a scalar product. ([m x 1]*[1 x m] = [1])
template<typename T, size_t m>
inline T operator*(const matrix<T,m,1>& lhs, const matrix<T,1,m>& rhs)
{
	T result(0);
	for(size_t i = 0; i < m; i++) {
		result += lhs(i,0) * rhs(0,i);
	}
	return result;
}

//! Returns the matrix scaled by the specified factor.
template<typename T, size_t m, size_t n>
inline matrix<T,m,n> operator*(const matrix<T,m,n>& mat, double factor)
{
	matrix<T,m,n> result(mat);
	result *= factor;
	return result;
}

//! Returns the matrix scaled by the specified factor.
template<typename T, size_t m, size_t n>
inline matrix<T,m,n> operator*(double factor, const matrix<T,m,n>& mat)
{
	matrix<T,m,n> result(mat);
	result *= factor;
	return result;
}

//! Returns the sum of the two matrices. Matrix dimensions must agree.
template<typename T, size_t m, size_t n>
inline matrix<T,m,n> operator+(const matrix<T,m,n>& lhs, const matrix<T,m,n>& rhs)
{
	matrix<T,m,n> result(lhs);
	result += rhs;
	return result;
}

//! Returns the difference of the two matrices. Matrix dimensions must agree.
template<typename T, size_t m, size_t n>
inline matrix<T,m,n> operator-(const matrix<T,m,n>& lhs, const matrix<T,m,n>& rhs)
{
	matrix<T,m,n> result(lhs);
	result -= rhs;
	return result;
}

#endif // MATRIX
