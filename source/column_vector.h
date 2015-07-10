/*
	linear_algebra_containers/column_vector header file
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

#ifndef COLUMN_VECTOR_H
#define COLUMN_VECTOR_H

#include <cstdlib>
#include <array>
#include <utility>
#include <cmath>
#include <algorithm>

#include "matrix.h"

/**
 * Column vector template with support of basic vector operations
 *
 * This template provides a basic column vector implementation for n-dimensional
 * vectors. It is based on the matrix template and similarily is not inteded to
 * be used in high performance computing with many dimensions. The class
 * feautres some convenience methods that extend the operations of the matrix
 * class with some vector operations.
 * @tparam T Type used for the entries of the matrix. Must support basic
 * aritmethic operations.
 * @tparam dim Number of rows of the vector.
 */
template<class T, size_t dim>
class column_vector : public matrix<T,dim,1>
{
public:
	//! The type of this vector
	typedef column_vector<T,dim> vector_type;

	/**
	 * @brief Construct an empty vector
	 *
	 * All entries are uninitialized if T is of fundamental type or constructed
	 * using their default constructor if T is a complex type.
	 */
	column_vector() = default;

	/**
	 * @brief Construct a vector with specified entries
	 *
	 * This constructor initializes the vector with the specified entries.
	 * It also allows using the braced initializer list syntax to initalize
	 * the vector, e.g.: column_vector<T,3> vec{{1,2,3}}.
	 */
	column_vector(const std::array<T,dim>& array)
		: matrix(array)
	{
	}

	/**
	 * @brief Construct a vector with specified entries
	 *
	 * This constructor initializes the vector with the specified entries.
	 * It also allows using the braced initializer list syntax to initalize
	 * the matrix, e.g.: column_vector<T,3> vec{{1,2,3}}.
	 * This constructor uses perfect forwarind to pass the an r-value to the
	 * constructor of the base class
	 */
	column_vector(std::array<T,dim>&& array)
		: matrix(std::forward<std::array<T,n*m>>(array))
	{
	}

	/**
	 * @brief Construct a vector from a matrix
	 *
	 * Constructs a column vector from a [m x 1] matrix by copying all values
	 * to the data array of the vector.
	 */
	column_vector(const matrix<T,dim,1>& mat)
		: matrix<T,dim,1>(mat)
	{
	}

	/**
	 * @brief Construct a vector from a matrix
	 *
	 * Constructs a column vector from a [m x 1] r-value matrix by forwarding
	 * all values to the constructor of the data array.
	 */
	column_vector(matrix<T,dim,1>&& mat)
		: matrix<T,dim,1>(std::forward<matrix<T,dim,1>>(mat))
	{
	}

	/**
	 * @brief Calculate euclidean norm
	 *
	 * Calculates the  length (the euclidean/2-norm) of the column vector.
	 * @return The length of the vector.
	 */
	T length() const
	{
		T result(0);
		for(size_t i = 0; i < dim; i++) {
			result += this->entries[i]*this->entries[i];
		}

		using std::sqrt;
		return sqrt(result);
	}

	/**
	 * @brief Calculate squared euclidean norm
	 *
	 * Calculates the length (the euclidean/2-norm) of the column vector without
	 * taking the square root of the inner product.
	 * @return The length of the vector squared.
	 */
	T lengthSquared() const
	{
		T result(0);
		for(size_t i = 0; i < dim; i++) {
			result += this->entries[i]*this->entries[i];
		}

		return result;
	}

	/**
	 * @brief Normalize the vector
	 *
	 * Normalizes the vector by dividing all entries of it by the length of the
	 * vector.
	 */
	void normalize()
	{
		(*this) *= (1/length());
	}

	/**
	 * @brief Create normalized copy
	 *
	 * Creates a normalized copy of the vector.
	 * @return The normalized copy of the vector.
	 */
	vector_type normalized() const
	{
		vector_type copy(*this);
		copy.normalize();
		return copy;
	}

	/**
	 * @brief Calculate the inner product
	 *
	 * This method calculates the inner product (dot product) of two column
	 * vectors.
	 * @param v1 The first vector.
	 * @param v2 The second vector.
	 * @return The inner product of the two vectors.
	 */
	static T dotProduct(const vector_type& v1, const vector_type& v2)
	{
		T result(0);
		for(size_t i = 0; i < dim; i++) {
			result += v1.entries[i] * v2.entries[i];
		}
		return result;
	}
};

#endif // COLUMN_VECTOR_H
