/*
	linear_algebra_containers/vector3 header file
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

#ifndef VECTOR3
#define VECTOR3

#include <utility>

#include "column_vector.h"

/**
 * Templated vector type for vectors in 3d space.
 *
 * @tparam T Type used for the entries of the vector.
 */
template<typename T>
class vector3 : public column_vector<T,3>
{

public:
	//! Constructs a 3d vector with uninitialized values.
	vector3() = default;

	//! Constructs a 3d vector with the specified values.
	vector3(T x, T y, T z)
	{
		this->entries[0] = x;
		this->entries[1] = y;
		this->entries[2] = z;
	}

	//! Constructs a 3d vector from the specified 3 entry vector.
	vector3(const column_vector<T,3>& vec)
		: column_vector<T,3>(vec)
	{
	}

	//! Constructs a new column vector from the specified 3 entry rvalue vector.
	vector3(column_vector<T,3>&& vec)
		: column_vector<T,3>(std::forward<column_vector<T,3>>(vec))
	{
	}

	//! Returns the x value of the vector.
	T x() const
	{
		return this->entries[0];
	}

	//! Returns the y value of the vector
	T y() const
	{
		return this->entries[1];
	}

	//! Returns the z value of the vector
	T z() const
	{
		return this->entries[2];
	}

	//! Sets the x value of the vector.
	void setX(T value)
	{
		this->entries[0] = value;
	}

	//! Sets the y value of the vector.
	void setY(T value)
	{
		this->entries[1] = value;
	}

	//! Sets the z value of the vector.
	void setZ(T value)
	{
		this->entries[2] = value;
	}

	//! Sets all values of the vector
	void set(T x, T y, T z)
	{
		this->entries[0] = x;
		this->entries[1] = y;
		this->entries[2] = z;
	}

	//! Calculates the cross product of two vectors.
	static vector3 crossProduct(const vector3& lhs, const vector3& rhs)
	{
		vector3 result;
		result[0] = lhs[1]*rhs[2] - lhs[2]*rhs[1];
		result[1] = lhs[2]*rhs[0] - lhs[0]*rhs[2];
		result[2] = lhs[0]*rhs[1] - lhs[1]*rhs[0];
		return result;
	}
};

#endif // VECTOR3
