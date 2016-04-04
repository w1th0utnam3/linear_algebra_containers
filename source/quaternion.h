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
#include <tuple>

#include "vector3.h"

/**
 * Quaternion template class with basic quaternion algebra operations
 *
 * This template provides a basic quaternion implementation. It represents
 * a standard quaternion: q = q0 + q1*i + q2*j + q3*k as a combination
 * of a scalar (q0) and a vector (q1, q2, q3).
 * Currently there are some operations implemented to help with rotation
 * operations of 3 dimensional vectors.
 *
 * @tparam T Type used for the entries of the quaternion.
 */
template<typename T>
class quaternion
{
protected:
	T _q0;				// Scalar component of the quaternion
	vector3<T> _qv;		// Vector component of the quaternion (q1, q2, q3)

public:
	// TODO: Different rotation orders
	// TODO: toMatrix
	// TODO: toEulerAngles

	//! Constructs an identity quaternion q(1 + 0*i + 0*j + 0*k).
	quaternion() : _q0(1), _qv(0,0,0) {}

	//! Constructs a quaternion q(q0 + q1*i + q2*j + q3*k).
	quaternion(T q0, T q1, T q2, T q3)
		: _q0(q0), _qv(q1,q2,q3) {}

	//! Constructs a quaternion q(q0, qv).
	quaternion(T q0, vector3<T> qv)
		: _q0(q0), _qv(qv) {}

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
		T q0q0 = _q0*_q0;
		if(q0q0 > 1) {
			axisOut->set(T(1), T(0), T(0));
			*angleOut = T(0);
		} else {
			using std::acos;
			using std::sqrt;
			*angleOut = 2*acos(_q0);
			T s = sqrt(1-q0q0);
			*axisOut = (1/s)*_qv;
		}
	}

	//! Returns a tuple containing the axis/angle representation of the quaternion. Quaternion must be normalized.
	std::tuple<vector3<T>,T> getAxisAndAngle()
	{
		std::tuple<vector3<T>,T> axisAngle;
		getAxisAndAngle(&std::get<vector3<T>>(axisAngle),&std::get<T>(axisAngle));
		return axisAngle;
	}

	//! Returns the conjugate of this quaternion.
	quaternion conjugated() const
	{
		return quaternion(_q0,-_qv);
	}

	//! Sets this quaternion to its conjugate.
	void conjugate()
	{
		_qv = -_qv;
	}

	//! Returns the dot product of two quaternions.
	static T dotProduct(const quaternion& p, const quaternion& q)
	{
		return p._q0*q._q0 + vector3<T>::dotProduct(p._qv,q._qv);
	}

	//! Returns the squared 2 norm of the quaternion.
	T normSquared() const
	{
		// return q0*q0 + qv.normSquared();
		return dotProduct(*this,*this);
	}

	//! Returns the 2 norm of the quaternion.
	T norm() const
	{
		using std::sqrt;
		return sqrt(this->normSquared());
	}

	//! Returns the normalized unit from this quaterion.
	quaternion normalized() const
	{
		T l = norm();
		return quaternion(_q0/l, (1/l)*_qv);
	}

	//! Normalizes this quaternion.
	void normalize()
	{
		double l = norm();
		_q0 = _q0/l;
		_qv *= (1/l);
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
		_q0 /= l2;
		_qv *= (1/l2);
	}

	//! Returns the scalar/real component of the quaternion
	T scalar() const
	{
		return _q0;
	}

	//! Returns the vector component of the quaternion
	vector3<T> vector() const
	{
		return _qv;
	}

	//! Returns the scalar component of the quaternion
	T q0() const
	{
		return _q0;
	}

	//! Returns the i coefficient of the quaternion
	T q1() const
	{
		return _qv.x();
	}

	//! Returns the j coefficient of the quaternion
	T q2() const
	{
		return _qv.y();
	}

	//! Returns the k coefficient of the quaternion
	T q3() const
	{
		return _qv.z();
	}

    //! Returns the transformation of the specified vector by this quaternion. Quaternion must be normalized!
    vector3<T> transform(const vector3<T>& v) const
    {
		return 2*_qv*(_qv.transposed()*v)-v*(_qv.transposed()*_qv)+(_q0*_q0)*v+2*_q0*(vector3<T>::crossProduct(_qv,v));
    }

	//! Returns the logarithm of the quaternion
	static quaternion log(const quaternion& p)
	{
		using std::log;
		using std::acos;

		auto lq = p.norm();
		auto lv = p._qv.norm();
		auto a = acos(p._q0 / lq)/lv;

		return quaternion(log(lq), a*p._qv);
	}

	//! Returns the exponential value of the quaternion
	static quaternion exp(const quaternion& p)
	{
		using std::exp;
		using std::cos;
		using std::sin;

		auto lv = p._qv.norm();
		auto s = sin(lv) / lv;
		return exp(p._q0)*quaternion(cos(lv), s*p._qv);
	}

	//! Returns the quaternion to the power of t
	static quaternion pow(const quaternion& p, T t)
	{
		return exp(t*log(p));
	}

    //! SO(3) group addition operator
	static quaternion composition(const quaternion& p, const quaternion& q)
    {
		return p*exp(0.5*q);
    }

    //! SO(3) group difference operator
	static quaternion difference(const quaternion& p, const quaternion& q)
    {
		return 2*log(q.inverse()*p);
    }

    //! Linear interpolation between two quaternions
	static quaternion slerp(const quaternion& p, const quaternion& q, T t)
    {
		return composition(p,t*(difference(q,p)));
    }

	//! Composition operator for quaternions
	friend quaternion operator*(const quaternion& p, const quaternion& q)
	{
		return quaternion(p._q0*q._q0 - p._qv.transposed()*q._qv, p._q0*q._qv + q._q0*p._qv + vector3<T>::crossProduct(p._qv, q._qv));
	}

	//! Scales a quaternion
	friend quaternion operator*(const T& n, const quaternion& p)
	{
		return quaternion(n*p._q0, n*p._qv);
	}

    //! Componentwise sum of two quaternions
	friend quaternion operator+(const quaternion& p, const quaternion& q)
	{
		return quaternion(p._q0 + q._q0, p._qv + q._qv);
	}

    //! Componentwise difference of two quaternions
	friend quaternion operator-(const quaternion& p, const quaternion& q)
	{
		return quaternion(p._q0 - q.q0, p._qv - q.qv);
	}

	//! Returns whether two quaternions have the same components
	friend bool operator==(const quaternion& lhs, const quaternion& rhs)
	{
		return ((lhs._q0 == rhs._q0) && (lhs._qv == rhs._qv));
	}

	//! Returns whether two quaternions do not have the same components
	friend bool operator!=(const quaternion& lhs, const quaternion& rhs)
	{
		return !(lhs == rhs);
	}
};

#endif // QUATERNION
