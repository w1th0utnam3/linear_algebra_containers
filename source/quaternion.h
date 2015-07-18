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
protected:
	T q0;				// Scalar component of the quaternion
	vector3<T> qv;		// Vector component of the quaternion

public:
	// TODO: Single member access
	// TODO: Different rotation orders
	// TODO: toMatrix
	// TODO: toEulerAngles

	//! Constructs an identity quaternion q(1 + 0*i + 0*j + 0*k).
	quaternion() : q0(1), qv(0,0,0) {}

	//! Constructs a quaternion q(q0 + q1*i + q2*j + q3*k).
	quaternion(T q0, T q1, T q2, T q3)
		: q0(q0), qv(q1,q2,q3) {}

	//! Constructs a quaternion q(q0, qv).
	quaternion(T q0, vector3<T> qv)
		: q0(q0), qv(qv) {}

	//! Creates a quaternion for the rotation about the specified angle around the specified axis. Axis must be normalized.
	static quaternion fromAxisAndAngle(const vector3<T>& axis, T angle)
	{
		using std::sin;
		using std::cos;
		T s = sin(angle/2);
		return quaternion(cos(angle/2), s*axis);
	}

	//! Creates a quaternion for the rotation about the specified angle around the specified axis. Axis must be normalized.
	static quaternion fromAxisAndAngle(T axisX, T axisY, T axisZ, T angle)
	{
		using std::sin;
		using std::cos;
		T s = sin(angle/2);
		return quaternion(cos(angle/2), axisX*s, axisY*s, axisZ*s);
	}

	//! Calculates the normalized axis/angle representation of the quaternion. Quaternion must be normalized.
	void getAxisAndAngle(vector3<T>* axisOut, T* angleOut)
	{
		T q0q0 = q0*q0;
		if(q0q0 > 1) {
			axisOut->set(T(1), T(0), T(0));
			angleOut = T(0);
		} else {
			using std::acos;
			using std::sqrt;
			angleOut = 2*acos(q0);
			T s = sqrt(1-q0q0);
			axisOut = (1/s)*qv;
		}
	}

	//! Returns the conjugate of this quaternion.
	quaternion conjugated() const
	{
		return quaternion(q0,-qv);
	}

	//! Sets this quaternion to its conjugate.
	void conjugate()
	{
		qv = -qv;
	}

	// Returns the 2 norm of the quaternion.
	T norm() const
	{
		using std::sqrt;
		return sqrt(q0*q0 + qv.normSquared());
	}

	// Returns the squared 2 norm of the quaternion.
	T normSquared() const
	{
		return q0*q0 + qv.normSquared();
	}

	//! Returns the normalized unit from this quaterion.
	quaternion normalized() const
	{
		T l = norm();
		return quaternion(q0/l, (1/l)*qv);
	}

	//! Normalizes this quaternion.
	void normalize()
	{
		double l = norm();
		q0 = q0/l;
		qv *= (1/l);
	}

	//! Returns the inverse of the quaternion
	quaternion inverse() const
	{
		return (1/this->normSquared())*(this->conjugated());
	}

	//! Inverts this quaternion
	quaternion invert()
	{
		conjugate();
		double l2 = normSquared();
		q0 /= l2;
		qv *= (1/l2);
	}

	//! Returns the scalar component of the quaternion
	T scalar() const
	{
		return q0;
	}

	//! Returns the vector component of the quaternion
	vector3<T> vector() const
	{
		return qv;
	}

	//! Returns the logarithm of the quaternion
	static quaternion log(const quaternion& q)
	{
		using std::log;
		using std::acos;

		auto lq = q.norm();
		auto lv = q.qv.norm();
		auto a = acos(q.q0 / lq)/lv;

		return quaternion(log(lq), a*q.qv);
	}

	//! Returns the exponential value of the quaternion
	static quaternion exp(const quaternion& q)
	{
		using std::exp;
		using std::cos;
		using std::sin;

		auto lv = q.qv.norm();
		auto s = sin(lv) / lv;
		return exp(q.q0)*quaternion(cos(lv), s*q.qv);
	}

	//! Returns the quaternion to the power of t
	static quaternion pow(const quaternion& q, T t)
	{
		return exp(t*log(q));
	}

	//! Linear interpolation between two quaternions
	static quaternion slerp(const quaternion& q1, const quaternion& q2, T t)
	{
		return q1+t*(q2-q1);
	}

	//! Composition operator for quaternions
	friend quaternion operator*(const quaternion& q, const quaternion& p)
	{
		return quaternion(q.q0*p.q0 - q.qv.transposed()*p.qv, q.q0*p.qv + p.q0*q.qv + vector3<T>::crossProduct(q.qv, p.qv));
	}

	//! Scales a quaternion
	friend quaternion operator*(const T& n, const quaternion& q)
	{
		return quaternion(n*q.q0, n*q.qv);
	}

	//! The SO(3) group addition operator (not componentwise sum!)
	friend quaternion operator+(const quaternion& q, const quaternion& v)
	{
		return q*exp(0.5*v);
	}

	//! The SO(3) group difference operator (not componentwise sum!)
	friend quaternion operator-(const quaternion& q, const quaternion& p)
	{
		return 2*log(p.inverse()*q);
	}

	//! Returns whether two quaternions have the same components
	friend bool operator==(const quaternion& lhs, const quaternion& rhs)
	{
		return ((lhs.q0 == rhs.q0) && (lhs.qv == rhs.qv));
	}

	//! Returns whether two quaternions do not have the same components
	friend bool operator!=(const quaternion& lhs, const quaternion& rhs)
	{
		return !(lhs == rhs);
	}
};

#endif // QUATERNION
