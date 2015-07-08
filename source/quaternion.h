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

#include "vector3.h"

/**
 * Template quaternion class with basic quaternion algebra operations.
 *
 * @tparam T Type used for the entries of the quaternion.
 */
template<typename T>
class quaternion
{
public:
	// TODO: Single member access
	// TODO: Different rotation orders
	// TODO: toMatrix
	// TODO: toEulerAngles

	//! Constructs an identity quaternion (1,0,0,0).
	quaternion() : w(1), x(0), y(0), z(0) {}

	//! Constructs a quaternion (w,x,y,z).
	quaternion(const T& w, const T& x, const T& y, const T& z)
		: w(w), x(x), y(y), z(z) {}

	//! Constructs a quaternion (w,x,y,z).
	quaternion(T&& w, T&& x, T&& y, T&& z)
		: w(std::move(w)), x(std::move(x)), y(std::move(y)), z(std::move(z)) {}

	//! Creates a quaternion for the rotation about the specified angle around the specified axis. Axis must be normalized.
	inline static quaternion fromAxisAndAngle(const vector3<T>& axis, const T& angle)
	{
		using std::sin;
		using std::cos;
		T s = sin(angle/2);
		return quaternion(cos(angle/2), axis.x()*s, axis.y()*s, axis.z()*s);
	}

	//! Creates a quaternion for the rotation about the specified angle around the specified axis. Axis must be normalized.
	inline static quaternion fromAxisAndAngle(const T& axisX, const T& axisY, const T& axisZ, const T& angle)
	{
		using std::sin;
		using std::cos;
		T s = sin(angle/2);
		return quaternion(cos(angle/2), axisX*s, axisY*s, axisZ*s);
	}

	//! Calculates the normalized axis/angle representation of the quaternion. Quaternion must be normalized.
	inline void getAxisAndAngle(vector3<T>* axisOut, T* angleOut)
	{
		T w2 = w*w;
		if(w2 > 1) {
			axisOut->setX(T(1));
			axisOut->setY(T(0));
			axisOut->setZ(T(0));
			angleOut = T(0);
		} else {
			using std::acos;
			using std::sqrt;
			angleOut = 2*acos(w);
			T s = sqrt(1-w2);
			axisOut->setX(x/s);
			axisOut->setY(y/s);
			axisOut->setZ(z/s);
		}
	}

	//! Returns the conjugate of this quaternion.
	inline quaternion conjugated() const
	{
		return quaternion(w,-x,-y,-z);
	}

	//! Sets this quaternion to its conjugate.
	inline void conjugate()
	{
		x = -x;
		y = -y;
		z = -z;
	}

	// Returns the length / norm of the quaternion.
	inline T length() const
	{
		using std::sqrt;
		return sqrt(w*w + x*x + y*y + z*z);
	}

	//! Returns the normalized unit from this quaterion.
	inline quaternion normalized() const
	{
		T l = length();
		return quaternion(w/l,x/l,y/l,z/l);
	}

	//! Normalizes this quaternion.
	inline void normalize()
	{
		double l = length();
		w = w/l;
		x = x/l;
		y = y/l;
		z = z/l;
	}

	//! Adds two quaternions.
	inline friend quaternion operator+(const quaternion& lhs, const quaternion& rhs)
	{
		return quaternion(lhs.w + rhs.w, lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
	}

	//! Subtracts two quaternions.
	inline friend quaternion operator-(const quaternion& lhs, const quaternion& rhs)
	{
		return quaternion(lhs.w - rhs.w, lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
	}

	//! Multiplies two quaternions.
	inline friend quaternion operator*(const quaternion& q1, const quaternion& q2)
	{
		return quaternion(q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z,
						  q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y,
						  q1.w*q2.y - q1.x*q2.z + q1.y*q2.w + q1.z*q2.x,
						  q1.w*q2.z + q1.x*q2.y - q1.y*q2.x + q1.z*q2.w);
	}

	//! Scales a quaternion.
	inline friend quaternion operator*(const T& factor, const quaternion& q)
	{
		return quaternion(q.w*factor, q.x*factor, q.y*factor, q.z*factor);
	}

	//! Scales a quaternion.
	inline friend const quaternion operator*(const quaternion& q, const T& factor)
	{
		return quaternion(q.w*factor, q.x*factor, q.y*factor, q.z*factor);
	}

public:
	T w;
	T x;
	T y;
	T z;
};

#endif // QUATERNION
