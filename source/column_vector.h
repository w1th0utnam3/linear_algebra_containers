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

namespace lin_algebra {

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
	 * @brief Construct a column vector
	 *
	 * Constructs a vector without initilization or a parameter pack
	 * for initialization depending on the arguments of the constructor.
	 */
	template<typename ...Ts>
	column_vector(Ts... values)
		: matrix{{values...}}
	{
	}

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
		return v1.transposed()*v2;
	}

	/**
	 * @brief Calculate squared euclidean norm
	 *
	 * Calculates the length (the euclidean/2-norm) of the column vector without
	 * taking the square root of the inner product.
	 * @return The length of the vector squared.
	 */
	T normSquared() const
	{
		return vector_type::dotProduct(*this,*this);
	}

	/**
	 * @brief Calculate squared euclidean norm
	 *
	 * Calculates the length (the euclidean/2-norm) of the specified column vector without
	 * taking the square root of the inner product.
	 * @param v The vector to claculate the squared length for.
	 * @return The length of the vector.
	 */
	static T normSquared(const vector_type& v) {
		return v.normSquared();
	}

	/**
	 * @brief Calculate euclidean norm
	 *
	 * Calculates the length (the euclidean/2-norm) of the column vector.
	 * @return The length of the vector.
	 */
	T norm() const
	{
		using std::sqrt;
		return sqrt(this->normSquared());
	}

	/**
	 * @brief Calculate euclidean norm
	 *
	 * Calculates the length (the euclidean/2-norm) of the specified column vector.
	 * @param v The vector to claculate the length for.
	 * @return The length of the vector.
	 */
	static T norm(const vector_type& v) {
		return v.norm();
	}

	/**
	 * @brief Normalize the vector
	 *
	 * Normalizes the vector by dividing all entries of it by the length of the
	 * vector.
	 */
	void normalize()
	{
		(*this) *= (1/norm());
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
	 * @brief Create normalized copy
	 *
	 * Creates a normalized copy of the specified vector.
	 * @param v The vector to create a normalized copy of.
	 * @return The normalized copy of the vector.
	 */
	static vector_type normalized(const vector_type& v)
	{
		return v.normalized();
	}
};

}

#endif // COLUMN_VECTOR_H
