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

#include <initializer_list>

#include "column_vector.h"

namespace lin_algebra {

/**
 * Vector template for 3d space
 *
 * This template provides some convenience methods on top of the
 * column_vector class which are appropriate for calculations in three
 * dimensions.
 * @tparam T Type used for the entries of the vector. Must support basic
 * aritmethic operations.
 */
template<typename T>
class Vector3 : public ColumnVector<T,3>
{

public:
	/**
	 * @brief Construct an empty vector
	 *
	 * All entries are uninitialized if T is of fundamental type or constructed
	 * using their default constructor if T is a complex type.
	 */
	Vector3() = default;

	//! Construct a vector with the specified values
	Vector3(T x, T y, T z) : Matrix{x,y,z} {}

	//! Implicit conversion from ColumnVector to Vector3
	Vector3(const ColumnVector<T,3>& v) : Matrix(v) {}

	/**
	 * @brief Return x value
	 * @return The x value of the vector
	 */
	T x() const
	{
		return this->entries_[0];
	}

	/**
	 * @brief Return y value
	 * @return The y value of the vector
	 */
	T y() const
	{
		return this->entries_[1];
	}

	/**
	 * @brief Return z value
	 * @return The z value of the vector
	 */
	T z() const
	{
		return this->entries_[2];
	}

	/**
	 * @brief Set x value
	 * @param value The value that the x-component should be set to.
	 */
	void setX(T value)
	{
		this->entries_[0] = value;
	}

	/**
	 * @brief Set y value
	 * @param value The value that the y-component should be set to.
	 */
	void setY(T value)
	{
		this->entries_[1] = value;
	}

	/**
	 * @brief Set z value
	 * @param value The value that the z-component should be set to.
	 */
	void setZ(T value)
	{
		this->entries_[2] = value;
	}

	/**
	 * @brief Set all values
	 * @param The value that the x-component should be set to.
	 * @param The value that the y-component should be set to.
	 * @param The value that the z-component should be set to.
	 */
	void set(T x, T y, T z)
	{
		this->entries_[0] = x;
		this->entries_[1] = y;
		this->entries_[2] = z;
	}

	/**
	 * @brief Calculate cross product
	 *
	 * Calculates the cross product of the two specified vectors.
	 * @param lhs The left vector.
	 * @param rhs The right vector.
	 * @return The vector normal to the two specified vectors.
	 */
	static Vector3 crossProduct(const Vector3& lhs, const Vector3& rhs)
	{
		Vector3 result;
		result[0] = lhs[1]*rhs[2] - lhs[2]*rhs[1];
		result[1] = lhs[2]*rhs[0] - lhs[0]*rhs[2];
		result[2] = lhs[0]*rhs[1] - lhs[1]*rhs[0];
		return result;
	}
};

}

#endif // VECTOR3
