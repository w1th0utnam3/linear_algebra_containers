/*
	linear_algebra_containers/quaternion header file
	MIT License

	Copyright (c) 2016 Fabian Löschner

	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE.
*/

#pragma once

#include "vector3.h"

#include <cmath>
#include <utility>
#include <tuple>

namespace lin_algebra {

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
class Quaternion
{
protected:
	T _q0;				// Scalar component of the quaternion
	Vector3<T> _qv;		// Vector component of the quaternion (q1, q2, q3)

public:
	// TODO: Different rotation orders
	// TODO: toMatrix
	// TODO: toEulerAngles

	//! Constructs an identity quaternion q(1 + 0*i + 0*j + 0*k).
	Quaternion() : _q0(1), _qv(0,0,0) {}

	//! Constructs a quaternion q(q0 + q1*i + q2*j + q3*k).
	Quaternion(T q0, T q1, T q2, T q3)
		: _q0(q0), _qv(q1,q2,q3) {}

	//! Constructs a quaternion q(q0, qv).
	Quaternion(T q0, Vector3<T> qv)
		: _q0(q0), _qv(qv) {}

	//! Creates a quaternion for the rotation about the specified angle around the specified axis. Axis must be normalized.
	static Quaternion fromAxisAndAngle(const Vector3<T>& axis, T angle)
	{
		using std::sin;
		using std::cos;
		T s = sin(angle/2);
		return Quaternion(cos(angle/2), s*axis);
	}

	//! Creates a quaternion for the rotation about the specified angle around the specified axis. Axis must be normalized.
	static Quaternion fromAxisAndAngle(T axisX, T axisY, T axisZ, T angle)
	{
		using std::sin;
		using std::cos;
		T s = sin(angle/2);
		return Quaternion(cos(angle/2), axisX*s, axisY*s, axisZ*s);
	}

	//! Calculates the normalized axis/angle representation of the quaternion. Quaternion must be normalized.
	void getAxisAndAngle(Vector3<T>* axisOut, T* angleOut)
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
	std::tuple<Vector3<T>,T> getAxisAndAngle()
	{
		std::tuple<Vector3<T>,T> axisAngle;
		getAxisAndAngle(&std::get<Vector3<T>>(axisAngle),&std::get<T>(axisAngle));
		return axisAngle;
	}

	//! Returns the conjugate of this quaternion.
	Quaternion conjugated() const
	{
		return Quaternion(_q0,-_qv);
	}

	//! Sets this quaternion to its conjugate.
	void conjugate()
	{
		_qv = -_qv;
	}

	//! Returns the dot product of two quaternions.
	static T dotProduct(const Quaternion& p, const Quaternion& q)
	{
		return p._q0*q._q0 + Vector3<T>::dotProduct(p._qv,q._qv);
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
	Quaternion normalized() const
	{
		T l = norm();
		return Quaternion(_q0/l, (1/l)*_qv);
	}

	//! Normalizes this quaternion.
	void normalize()
	{
		double l = norm();
		_q0 = _q0/l;
		_qv *= (1/l);
	}

	//! Returns the inverse of the quaternion
	Quaternion inverse() const
	{
		return (1/this->normSquared())*(this->conjugated());
	}

	//! Inverts this quaternion
	Quaternion invert()
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
	Vector3<T> vector() const
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
	Vector3<T> transform(const Vector3<T>& v) const
	{
		return 2*_qv*(_qv.transposed()*v)-v*(_qv.normSquared())+(_q0*_q0)*v+2*_q0*(Vector3<T>::crossProduct(_qv,v));
	}

	//! Returns the logarithm of the quaternion
	static Quaternion log(const Quaternion& p)
	{
		using std::log;
		using std::acos;

		auto lq = p.norm();
		auto lv = p._qv.norm();
		auto a = acos(p._q0 / lq)/lv;

		return Quaternion(log(lq), a*p._qv);
	}

	//! Returns the exponential value of the quaternion
	static Quaternion exp(const Quaternion& p)
	{
		using std::exp;
		using std::cos;
		using std::sin;

		auto lv = p._qv.norm();
		auto s = sin(lv) / lv;
		return exp(p._q0)*Quaternion(cos(lv), s*p._qv);
	}

	//! Returns the quaternion to the power of t
	static Quaternion pow(const Quaternion& p, T t)
	{
		return exp(t*log(p));
	}

	//! SO(3) group addition operator
	static Quaternion composition(const Quaternion& p, const Quaternion& q)
	{
		return p*exp(0.5*q);
	}

	//! SO(3) group difference operator
	static Quaternion difference(const Quaternion& p, const Quaternion& q)
	{
		return 2*log(q.inverse()*p);
	}

	//! Linear interpolation between two quaternions
	static Quaternion slerp(const Quaternion& p, const Quaternion& q, T t)
	{
		return composition(p,t*(difference(q,p)));
	}

	//! Composition operator for quaternions
	friend Quaternion operator*(const Quaternion& p, const Quaternion& q)
	{
		return Quaternion(p._q0*q._q0 - p._qv.transposed()*q._qv, p._q0*q._qv + q._q0*p._qv + Vector3<T>::crossProduct(p._qv, q._qv));
	}

	//! Scales a quaternion
	friend Quaternion operator*(const T& n, const Quaternion& p)
	{
		return Quaternion(n*p._q0, n*p._qv);
	}

	//! Componentwise sum of two quaternions
	friend Quaternion operator+(const Quaternion& p, const Quaternion& q)
	{
		return Quaternion(p._q0 + q._q0, p._qv + q._qv);
	}

	//! Componentwise difference of two quaternions
	friend Quaternion operator-(const Quaternion& p, const Quaternion& q)
	{
		return Quaternion(p._q0 - q.q0, p._qv - q.qv);
	}

	//! Returns whether two quaternions have the same components
	friend bool operator==(const Quaternion& lhs, const Quaternion& rhs)
	{
		return ((lhs._q0 == rhs._q0) && (lhs._qv == rhs._qv));
	}

	//! Returns whether two quaternions do not have the same components
	friend bool operator!=(const Quaternion& lhs, const Quaternion& rhs)
	{
		return !(lhs == rhs);
	}
};

//! Prints the quaternion to the specified stream
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const Quaternion<T>& quat)
{
	os << "[";
	os << quat.q0() << ";" << quat.q1() << ";" << quat.q2() << ";" << quat.q3() << ";";
	os << "]";

	return os;
}

}
