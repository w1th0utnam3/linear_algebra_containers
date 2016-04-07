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
class vector3 : public column_vector<T,3>
{

public:
	/**
	 * @brief Construct an empty vector
	 *
	 * All entries are uninitialized if T is of fundamental type or constructed
	 * using their default constructor if T is a complex type.
	 */
	vector3() = default;

	/**
	 * @brief Construct a vector with values
	 *
	 * Constructs a vector initialized to the specified values.
	 */
	vector3(T x, T y, T z)
		: column_vector<T,3>{x,y,z}
	{
	}

	/**
	* @brief Construct a vector3 from column vector
	*
	* Constructs a vector3 from the specified 3 dimensional column vector.
	*/
	vector3(const column_vector<T,3>& vec)
		: column_vector<T,3>(vec)
	{
	}

	/**
	* @brief Construct a vector3 from matrix
	*
	* Constructs a vector3 from the specified 3x1 dimensional matrix.
	*/
	vector3(const matrix<T,3,1>& mat)
		: column_vector<T,3>(mat)
	{
	}

	/**
	* @brief Construct a vector3 from column vector
	*
	* Constructs a vector3 from the specified 3 dimensional r-value column
	* vector.
	*/
	vector3(column_vector<T,3>&& vec)
		: column_vector<T,3>(std::forward<column_vector<T,3>>(vec))
	{
	}

	/**
	 * @brief Return x value
	 * @return The x value of the vector
	 */
	T x() const
	{
		return this->_entries[0];
	}

	/**
	 * @brief Return y value
	 * @return The y value of the vector
	 */
	T y() const
	{
		return this->_entries[1];
	}

	/**
	 * @brief Return z value
	 * @return The z value of the vector
	 */
	T z() const
	{
		return this->_entries[2];
	}

	/**
	 * @brief Set x value
	 * @param value The value that the x-component should be set to.
	 */
	void setX(T value)
	{
		this->_entries[0] = value;
	}

	/**
	 * @brief Set y value
	 * @param value The value that the y-component should be set to.
	 */
	void setY(T value)
	{
		this->_entries[1] = value;
	}

	/**
	 * @brief Set z value
	 * @param value The value that the z-component should be set to.
	 */
	void setZ(T value)
	{
		this->_entries[2] = value;
	}

	/**
	 * @brief Set all values
	 * @param The value that the x-component should be set to.
	 * @param The value that the y-component should be set to.
	 * @param The value that the z-component should be set to.
	 */
	void set(T x, T y, T z)
	{
		this->_entries[0] = x;
		this->_entries[1] = y;
		this->_entries[2] = z;
	}

	/**
	 * @brief Calculate cross product
	 *
	 * Calculates the cross product of the two specified vectors.
	 * @param lhs The left vector.
	 * @param rhs The right vector.
	 * @return The vector normal to the two specified vectors.
	 */
	static vector3 crossProduct(const vector3& lhs, const vector3& rhs)
	{
		vector3 result;
		result[0] = lhs[1]*rhs[2] - lhs[2]*rhs[1];
		result[1] = lhs[2]*rhs[0] - lhs[0]*rhs[2];
		result[2] = lhs[0]*rhs[1] - lhs[1]*rhs[0];
		return result;
	}
};

}

#endif // VECTOR3
