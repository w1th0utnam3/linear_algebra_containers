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

#include "matrixbase.h"
#include "matrix.h"

namespace lin_algebra {

/**
 * Column vector template with support of basic vector operations
 *
 * This template provides a basic column vector implementation for n-dimensional
 * vectors. It is a partial specialization of the matrix template and similarily
 * is not inteded to be used in high performance computing with many dimensions.
 * The class feautres some convenience methods that extend the operations of the
 * matrix class with some vector operations.
 * @tparam T Type used for the entries of the vector. Must support basic
 * aritmethic operations.
 * @tparam dim Number of rows of the vector.
 */
template<class T, size_t dim>
class Matrix<T,dim,1> : public MatrixBase<T,dim,1>
{

public:
	//! Type of this vector
	typedef Matrix<T,dim,1> vector_type;
	//! Type of the transposed vector
	typedef Matrix<T,1,dim> transposed_vector_type;

	/**
	 * @brief Construct a column vector
	 *
	 * Constructs a vector either without initilization or a parameter pack
	 * for initialization depending on the arguments of the constructor.
	 */
	template<typename ...Ts>
	Matrix(Ts... values)
		: MatrixBase{values...}
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
		T product{T(0)};
		for(size_t i = 0; i < rows; i++) {
			product += v1[i]*v2[i];
		}
		return product;
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
		T norm{T(0)};
		for(size_t i = 0; i < rows; i++) {
			norm += this->entries_[i]*this->entries_[i];
		}
		return norm;
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

	//! Returns the transposed vector
	transposed_vector_type transposed() const
	{
		return transposed_vector_type{this->entries_};
	}

	//! Adds the right vector to the left.
	vector_type operator+=(const vector_type& rhs)
	{
		for(size_t i = 0; i < rows; i++) {
			this->entries_[i] += rhs.entries_[i];
		}
		return *this;
	}

	//! Substracts the right vector from the left.
	vector_type operator-=(const vector_type& rhs)
	{
		for(size_t i = 0; i < rows; i++) {
			this->entries_[i] -= rhs.entries_[i];
		}
		return *this;
	}

	//! Scales the vector by the specified factor.
	vector_type operator*=(double factor)
	{
		for(size_t i = 0; i < rows; i++) {
			this->entries_[i] *= factor;
		}
		return *this;
	}

	//! Returns the vector scaled by the specified factor.
	friend vector_type operator*(const vector_type& mat, double factor)
	{
		vector_type result(mat);
		result *= factor;
		return result;
	}

	//! Returns the vector scaled by the specified factor.
	friend vector_type operator*(double factor, const vector_type& mat)
	{
		vector_type result(mat);
		result *= factor;
		return result;
	}

	//! Returns the sum of the two matrices. Vector dimensions must agree.
	friend vector_type operator+(const vector_type& lhs, const vector_type& rhs)
	{
		vector_type result(lhs);
		result += rhs;
		return result;
	}

	//! Returns the difference of the two matrices. Vector dimensions must agree.
	friend vector_type operator-(const vector_type& lhs, const vector_type& rhs)
	{
		vector_type result(lhs);
		result -= rhs;
		return result;
	}

	//! Returns the negated vector.
	friend vector_type operator-(const vector_type& in)
	{
		return T(-1)*vector_type(in);
	}
};

//! Column vector template alias
template<typename T, size_t dim>
using ColumnVector = Matrix<T,dim,1>;

}

#endif // COLUMN_VECTOR_H
