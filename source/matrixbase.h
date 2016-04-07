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

#ifndef MATRIXBASE_H
#define MATRIXBASE_H

#include <array>

namespace lin_algebra {

/**
* Base class for matrices and vectors
*
* This template provides a basic matrix implementation for m x n matrices. Storage
* of the entries is column-major. For further explenation see Matrix subclass.
* This class specifies only methods that don't return matrix types so that
* subclassing and partial specialization allows to return specialized types
* (i.e. vec+vec=vec, mat*vec=vec to allow (mat*vec).norm()=scalar etc.).
*/
template<typename T, size_t row_count_param, size_t column_count_param>
class MatrixBase
{
public:
	//! Number of rows of this matrix type
	static constexpr size_t rows = row_count_param;
	//! Number of columns of this matrix type
	static constexpr size_t cols = column_count_param;

protected:
	//! Array storing the matrix entries
	std::array<T,rows*cols> entries_;

public:
	//! Constructs a matrix (either unintialized if called without arguments or initialized with the specified values)
	template<typename ...Ts>
	MatrixBase(Ts... values) : entries_{values...} {}

	//! Returns the element index of the specified coordinates
	static constexpr size_t index(size_t row, size_t column) { return row + column*rows; }

	//! Returns a reference to the entry at the specified coordinates
	T& operator()(size_t row, size_t column) { return entries_[MatrixBase::index(row,column)]; }
	//! Returns a const reference to the entry at the specified coordinates
	const T& operator()(size_t row, size_t column) const { return entries_[MatrixBase::index(row,column)]; }

	//! Returns a reference to the i-th element stored in the matrix (column major)
	T& operator[](size_t i) { return entries_[i]; }
	//! Returns a const-reference to the i-th element stored in the matrix (column major)
	const T& operator[](size_t i) const { return entries_[i]; }

	//! Returns a pointer to the underlying array (column major)
	T* data() { return entries_.data(); }
	//! Returns a const-pointer to the underlying array (column major)
	const T* data() const { return entries_.data(); }

	//! Sets all entries to the specified value
	void fill(const T& val) { for(auto& v : entries_) v = val; }
	//! Sets all entries to zero
	void zeros() { this->fill(T(0)); }

	//! Compares the matrices elementwise for equality
	friend bool operator==(const MatrixBase& lhs, const MatrixBase& rhs) { return lhs.entries_ == rhs.entries_; }
	//! Compares the matrices elementwise for inequality
	friend bool operator!=(const MatrixBase& lhs, const MatrixBase& rhs) { return !(lhs == rhs); }
};

//! Prints the matrix to the specified stream
template<typename T, size_t m, size_t n>
inline std::ostream& operator<<(std::ostream& os, const MatrixBase<T,m,n>& mat)
{
	os << "[";
	for(size_t i = 0; i < m; i++) {
		for(size_t j = 0; j < n-1; j++) os << mat(i,j) << " ";
		os << mat(i,n-1) << ";";
		if(i < m-1) os << " ";
	}
	os << "]";

	return os;
}

}

#endif // MATRIXBASE_H
