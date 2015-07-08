/*
	linear_algebra_containers/quaternion header file
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

#ifndef QUATERNION
#define QUATERNION

#include <cmath>
#include <utility>

/**
 * Template quaternion class with basic quaternion algebra operations.
 *
 * @tparam T Type used for the entries of the quaternion.
 */
template<typename T>
class Quaternion
{
public:
	// TODO: Single member access
	// TODO: fromAxisAngle
	// TODO: Different rotation orders
	// TODO: toMatrix
	// TODO: toEulerAngles
	// TODO: toAxisAngle

	//! Constructs an identity quaternion (1,0,0,0).
	Quaternion() : w(1), x(0), y(0), z(0) {}

	//! Constructs a quaternion (w,x,y,z).
	Quaternion(const T& w, const T& x, const T& y, const T& z)
		: w(w), x(x), y(y), z(z) {}

	//! Constructs a quaternion (w,x,y,z).
	Quaternion(T&& w, T&& x, T&& y, T&& z)
		: w(std::move(w)), x(std::move(x)), y(std::move(y)), z(std::move(z)) {}

	//! Returns the conjugate of this quaternion
	Quaternion conjugated() const
	{
		return Quaternion(w,-x,-y,-z);
	}

	//! Sets this quaternion to its conjugate
	void conjugate()
	{
		x = -x;
		y = -y;
		z = -z;
	}

	// Returns the length / norm of the quaternion
	T length() const
	{
		using std::sqrt;
		return sqrt(w*w + x*x + y*y + z*z);
	}

	//! Returns the normalized unit from this quaterion
	Quaternion normalized() const
	{
		T l = length();
		return Quaternion(w/l,x/l,y/l,z/l);
	}

	//! Normalizes this quaternion
	void normalize()
	{
		double l = length();
		w = w/l;
		x = x/l;
		y = y/l;
		z = z/l;
	}

	friend Quaternion operator+(const Quaternion& lhs, const Quaternion& rhs)
	{
		return Quaternion(lhs.w + rhs.w, lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
	}

	friend Quaternion operator-(const Quaternion& lhs, const Quaternion& rhs)
	{
		return Quaternion(lhs.w - rhs.w, lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
	}

	friend Quaternion operator*(const Quaternion& q1, const Quaternion& q2)
	{
		return Quaternion(q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z,
						  q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y,
						  q1.w*q2.y - q1.x*q2.z + q1.y*q2.w + q1.z*q2.x,
						  q1.w*q2.z + q1.x*q2.y - q1.y*q2.x + q1.z*q2.w);
	}

	friend Quaternion operator*(const T& factor, const Quaternion& q)
	{
		return Quaternion(q.w*factor, q.x*factor, q.y*factor, q.z*factor);
	}

	friend const Quaternion operator*(const Quaternion& q, const T& factor)
	{
		return Quaternion(q.w*factor, q.x*factor, q.y*factor, q.z*factor);
	}

public:
	T w;
	T x;
	T y;
	T z;
};

#endif // QUATERNION
