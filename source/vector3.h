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

#include <cmath>

#include "matrixbase.h"
#include "matrix.h"

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
class Matrix<T,3,1> : public MatrixBase<T,3,1>
{
private:
	typedef MatrixBase<T,3,1> MatrixBaseType;

public:
	//! Type of this vector
	typedef Matrix<T,3,1> VectorType;
	//! Type of the transposed vector
	typedef Matrix<T,1,3> TransposedVectorType;

	/**
	 * @brief Construct an empty vector
	 *
	 * All entries are uninitialized if T is of fundamental type or constructed
	 * using their default constructor if T is a complex type.
	 */
	Matrix() = default;

	//! Construct a vector with the specified values
	Matrix(T x, T y, T z)
		: MatrixBaseType(x,y,z)
	{
	}

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
	static VectorType crossProduct(const VectorType& lhs, const VectorType& rhs)
	{
		VectorType result;
		result[0] = lhs[1]*rhs[2] - lhs[2]*rhs[1];
		result[1] = lhs[2]*rhs[0] - lhs[0]*rhs[2];
		result[2] = lhs[0]*rhs[1] - lhs[1]*rhs[0];
		return result;
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
	static T dotProduct(const VectorType& v1, const VectorType& v2)
	{
		return (  v1[0]*v2[0]
				+ v1[1]*v2[1]
				+ v1[2]*v2[2] );
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
		return (  this->entries_[0]*this->entries_[0]
				+ this->entries_[1]*this->entries_[1]
				+ this->entries_[2]*this->entries_[2] );
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
	 * @brief Normalize the vector
	 *
	 * Normalizes the vector by dividing all entries of it by the length of the
	 * vector.
	 */
	VectorType& normalize()
	{
		(*this) *= (1/norm());
		return *this;
	}

	/**
	 * @brief Create normalized copy
	 *
	 * Creates a normalized copy of the vector.
	 * @return The normalized copy of the vector.
	 */
	VectorType normalized() const
	{
		VectorType copy(*this);
		copy.normalize();
		return copy;
	}

	//! Returns the transposed vector
	TransposedVectorType transposed() const
	{
		return TransposedVectorType{this->entries_};
	}

	//! Adds the right vector to the left.
	VectorType operator+=(const VectorType& rhs)
	{
		this->entries_[0] += rhs.entries_[0];
		this->entries_[1] += rhs.entries_[1];
		this->entries_[2] += rhs.entries_[2];
		return *this;
	}

	//! Substracts the right vector from the left.
	VectorType operator-=(const VectorType& rhs)
	{
		this->entries_[0] -= rhs.entries_[0];
		this->entries_[1] -= rhs.entries_[1];
		this->entries_[2] -= rhs.entries_[2];
		return *this;
	}

	//! Scales the vector by the specified factor.
	VectorType operator*=(double factor)
	{
		this->entries_[0] *= factor;
		this->entries_[1] *= factor;
		this->entries_[2] *= factor;
		return *this;
	}

	//! Returns the vector scaled by the specified factor.
	friend VectorType operator*(const VectorType& mat, double factor)
	{
		VectorType result(mat);
		result *= factor;
		return result;
	}

	//! Returns the vector scaled by the specified factor.
	friend VectorType operator*(double factor, const VectorType& mat)
	{
		VectorType result(mat);
		result *= factor;
		return result;
	}

	//! Returns the sum of the two matrices. Vector dimensions must agree.
	friend VectorType operator+(const VectorType& lhs, const VectorType& rhs)
	{
		VectorType result(lhs);
		result += rhs;
		return result;
	}

	//! Returns the difference of the two matrices. Vector dimensions must agree.
	friend VectorType operator-(const VectorType& lhs, const VectorType& rhs)
	{
		VectorType result(lhs);
		result -= rhs;
		return result;
	}

	//! Returns the negated vector.
	friend VectorType operator-(const VectorType& in)
	{
		return T(-1)*VectorType(in);
	}
};

//! Vector template alias
template<typename T>
using Vector3 = Matrix<T,3,1>;

}

#endif // VECTOR3
