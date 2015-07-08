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
 * Template column vector type with support of basic vector operations.
 *
 * @tparam T Type used for the entries of the vector.
 * @tparam dim Number of rows of the vector.
 */
template<class T, size_t dim>
class column_vector : public matrix<T,dim,1>
{
public:
	//! The type of this vector
	typedef column_vector<T,dim> vector_type;

	//! Constructs a vector with uninitialized values.
	column_vector() = default;

	//! Constructs a new column vector from the specified [m x 1] matrix
	inline column_vector(const matrix<T,dim,1>& mat)
		: matrix<T,dim,1>(mat)
	{
	}

	//! Constructs a new column vector from the specified [m x 1] rvalue matrix
	inline column_vector(matrix<T,dim,1>&& mat)
		: matrix<T,dim,1>(std::forward<matrix<T,dim,1>>(mat))
	{
	}

	//! Returns the length (2 norm) of the vector.
	inline T length() const
	{
		T result(0);
		for(size_t i = 0; i < dim; i++) {
			result += this->entries[i]*this->entries[i];
		}

		using std::sqrt;
		return sqrt(result);
	}

	//! Returns the squared 2 norm of the vector.
	inline T lengthSquared() const
	{
		T result(0);
		for(size_t i = 0; i < dim; i++) {
			result += this->entries[i]*this->entries[i];
		}

		return result;
	}

	//! Normalizes the vector by dividing all elements by the vectors length.
	inline void normalize()
	{
		(*this) *= (1/length());
	}

	//! Returns a normalized copy of this vector.
	inline vector_type normalized() const
	{
		vector_type copy(*this);
		copy.normalize();
		return copy;
	}

	//! Calculates the dot product of two column vectors. Vector dimensions have to agree.
	inline static T dotProduct(const vector_type& v1, const vector_type& v2)
	{
		T result(0);
		for(size_t i = 0; i < dim; i++) {
			result += v1.entries[i] * v2.entries[i];
		}
		return result;
	}
};

#endif // COLUMN_VECTOR_H
