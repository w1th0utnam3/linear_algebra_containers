//	MIT License
//
//	Copyright (c) 2016 Fabian Löschner
//
//	Permission is hereby granted, free of charge, to any person obtaining a copy
//	of this software and associated documentation files (the "Software"), to deal
//	in the Software without restriction, including without limitation the rights
//	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//	copies of the Software, and to permit persons to whom the Software is
//	furnished to do so, subject to the following conditions:
//
//	The above copyright notice and this permission notice shall be included in all
//	copies or substantial portions of the Software.
//
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//	SOFTWARE.

#include "catch.hpp"

#include "matrix.h"
#include "column_vector.h"
#include "vector3.h"
#include "quaternion.h"

using namespace lin_algebra;

template<typename T>
static T relativeError(T value, T approxValue)
{
	using std::abs;
	return abs(T(1) - (value / approxValue));
}

TEST_CASE("Testing Matrix")
{
	// TODO: Test data
	// TODO: Test const
	// TODO: Test transposed
	// TODO: Test comparison operators
	// TODO: Test negation

	typedef Matrix<double, 4, 4> mat4x4d;
	typedef Matrix<double, 4, 2> mat4x2d;
	typedef Matrix<double, 2, 2> mat2x2d;

	mat4x4d mat{ 1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16. };

	SECTION("Testing dimension size")
	{
		REQUIRE(mat4x4d::rows == 4);
		REQUIRE(mat4x4d::cols == 4);
	}

	SECTION("Testing subscript operator")
	{
		for (int i = 0; i < 16; i++) {
			REQUIRE(mat[i] == i + 1);
		}
	}

	SECTION("Testing ()-operator")
	{
		REQUIRE(mat(0, 0) == 1);
		REQUIRE(mat(1, 0) == 2);
		REQUIRE(mat(2, 0) == 3);
		REQUIRE(mat(3, 0) == 4);
		REQUIRE(mat(0, 1) == 5);
		REQUIRE(mat(1, 1) == 6);
		REQUIRE(mat(2, 1) == 7);
		REQUIRE(mat(3, 1) == 8);
		REQUIRE(mat(0, 2) == 9);
		REQUIRE(mat(1, 2) == 10);
		REQUIRE(mat(2, 2) == 11);
		REQUIRE(mat(3, 2) == 12);
		REQUIRE(mat(0, 3) == 13);
		REQUIRE(mat(1, 3) == 14);
		REQUIRE(mat(2, 3) == 15);
		REQUIRE(mat(3, 3) == 16);
	}

	SECTION("Testing fill()")
	{
		const double value = 12.345;
		mat.fill(value);

		for (int i = 0; i < 16; i++) {
			REQUIRE(mat[i] == value);
		}
	}

	SECTION("Testing zeros()")
	{
		mat.zeros();

		for (int i = 0; i < 16; i++) {
			REQUIRE(mat[i] == 0);
		}
	}

	SECTION("Testing toIdentity()")
	{
		mat.toIdentity();
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				if (i == j) {
					REQUIRE(mat(i, j) == 1);
				} else {
					REQUIRE(mat(i, j) == 0);
				}
			}
		}
	}

	SECTION("Testing matrix multiplication")
	{
		for (int i = 0; i < 16; i++) {
			mat[i] = i + 1;
		}

		mat4x2d mat2;
		for (int i = 0; i < 8; i++) {
			mat2[i] = i + 1;
		}

		const auto result = mat*mat2;

		REQUIRE((std::is_same<decltype(result), const mat4x2d>::value == true));
		REQUIRE(result(0, 0) == 90);
		REQUIRE(result(1, 0) == 100);
		REQUIRE(result(2, 0) == 110);
		REQUIRE(result(3, 0) == 120);
		REQUIRE(result(0, 1) == 202);
		REQUIRE(result(1, 1) == 228);
		REQUIRE(result(2, 1) == 254);
		REQUIRE(result(3, 1) == 280);
	}

	SECTION("Testing matrix scaling")
	{
		Matrix<double, 2, 2> mat2;
		mat2.fill(3);
		mat2 *= 2;

		REQUIRE(mat2[0] == 6);
		REQUIRE(mat2[1] == 6);
		REQUIRE(mat2[2] == 6);
		REQUIRE(mat2[3] == 6);
	}

	SECTION("Testing matrix sum/difference")
	{
		mat2x2d mat1, mat2;
		mat1.fill(2);
		mat2.fill(6);

		mat2 -= mat1;
		REQUIRE(mat2[0] == 4);
		REQUIRE(mat2[1] == 4);
		REQUIRE(mat2[2] == 4);
		REQUIRE(mat2[3] == 4);

		mat1.fill(3);
		mat1 += mat2;
		REQUIRE(mat1[0] == 7);
		REQUIRE(mat1[1] == 7);
		REQUIRE(mat1[2] == 7);
		REQUIRE(mat1[3] == 7);

		mat2.toIdentity();
		const auto mat6 = mat1 - mat2;
		REQUIRE(mat6(0, 0) == 6);
		REQUIRE(mat6(1, 0) == 7);
		REQUIRE(mat6(0, 1) == 7);
		REQUIRE(mat6(1, 1) == 6);
	}

	SECTION("Testing matrix comparison")
	{
		mat4x4d mat2;
		mat2.toIdentity();

		REQUIRE(mat != mat2);
		REQUIRE(!(mat2 == mat));
	}

	SECTION("Testing matrix copy/move constructors")
	{
		mat4x4d mat2(mat);
		REQUIRE(mat2 == mat);

		mat.fill(0);
		mat4x4d mat3(mat4x4d().fill(0));
		REQUIRE(mat3 == mat);

		mat.toIdentity();
		mat4x4d mat4(std::move(mat4x4d().toIdentity()));
		REQUIRE(mat4 == mat);
	}
}

TEST_CASE("Testing ColumnVector")
{
	typedef ColumnVector<double, 4> vec4d;
	vec4d v1, v2;

	SECTION("Testing initializer list")
	{
		vec4d ref;
		ref.fill(4.5);

		vec4d test1(4.5, 4.5, 4.5, 4.5);
		vec4d test2 = { 4.5,4.5,4.5,4.5 };

		REQUIRE(test1 == ref);
		REQUIRE(test2 == ref);
	}

	SECTION("Testing fill() and array subscript operator")
	{
		v1.fill(2.5);
		REQUIRE(v1[0] == 2.5);
		REQUIRE(v1[1] == 2.5);
		REQUIRE(v1[2] == 2.5);
		REQUIRE(v1[3] == 2.5);
	}

	SECTION("Testing zeors()")
	{
		v2.zeros();
		REQUIRE(v2[0] == 0);
		REQUIRE(v2[1] == 0);
		REQUIRE(v2[2] == 0);
		REQUIRE(v2[3] == 0);
	}

	SECTION("Testing norm()")
	{
		v1.fill(3);
		REQUIRE(v1.norm() == sqrt(36));
	}

	SECTION("Testing normalize()")
	{
		v1.fill(3);
		v1.normalize();
		REQUIRE(v1.norm() == 1);
	}

	SECTION("Testing normalized()")
	{
		v1.fill(3);
		REQUIRE(v1.normalized().norm() == 1);
	}

	SECTION("Testing dotProduct()")
	{
		v1[0] = 1;
		v1[1] = 2;
		v1[2] = 3;
		v1[3] = 4;
		v2[0] = 5;
		v2[1] = 6;
		v2[2] = 7;
		v2[3] = 8;
		REQUIRE(vec4d::dotProduct(v1, v2) == 70);
	}

	SECTION("Testing normSquared()")
	{
		v1.fill(2);
		REQUIRE(v1.normSquared() == 16);
	}

	SECTION("Testing implicit conversion between matrix and column vector")
	{
		v1.fill(2); v2.fill(3);
		REQUIRE(vec4d::dotProduct((v1.transposed().transposed()), v2) == 24);

		auto t = v1.transposed().transposed();
		REQUIRE(vec4d::dotProduct(t, v2) == 24);

		Matrix<double, 4, 4> mat = Matrix<double, 4, 4>::createIdentity();
		REQUIRE(vec4d::dotProduct(mat*t, v2) == 24);
	}
}

TEST_CASE("Testing Vector3")
{
	typedef Vector3<double> vec3d;
	typedef Matrix<double, 3, 3> mat3x3d;

	vec3d v1, v2;

	SECTION("Testing constructor")
	{
		vec3d test1(0.1, 312.112, 77);
		REQUIRE(test1[0] == 0.1);
		REQUIRE(test1[1] == 312.112);
		REQUIRE(test1[2] == 77);

		const double k = 77;
		vec3d test2(0.1, 312.112, k);
		REQUIRE(test2[0] == 0.1);
		REQUIRE(test2[1] == 312.112);
		REQUIRE(test2[2] == k);
	}

	SECTION("Testing initializer list")
	{
		v1.fill(4);
		v2 = { 4,4,4 };
		REQUIRE(v1 == v2);
	}

	SECTION("Testing getters")
	{
		v1[0] = 3.3;
		v1[1] = 4.4;
		v1[2] = 5.5;
		REQUIRE(v1.x() == 3.3);
		REQUIRE(v1.y() == 4.4);
		REQUIRE(v1.z() == 5.5);
	}

	SECTION("Testing setters")
	{
		const double d = 22.2;
		v2.setX(1);
		v2.setY(d);
		v2.setZ(22);
		REQUIRE(v2.x() == 1);
		REQUIRE(v2.y() == d);
		REQUIRE(v2.z() == 22);
	}

	SECTION("Testing crossProduct()")
	{
		v1.setX(1);
		v1.setY(2);
		v1.setZ(3);
		v2.setX(3);
		v2.setY(4);
		v2.setZ(5);
		auto r = vec3d::crossProduct(v1, v2);
		REQUIRE(r.x() == -2);
		REQUIRE(r.y() == 4);
		REQUIRE(r.z() == -2);
	}

	SECTION("Testing conversions")
	{
		mat3x3d mat = mat3x3d::createIdentity();
		v1 = { 1,2,3 };
		REQUIRE((mat*v1).z() == 3);
	}
}

TEST_CASE("Testing Quaternion")
{
	typedef Quaternion<double> quatd;
	typedef Vector3<double> vec3d;

	quatd q;

	SECTION("Testing constructor")
	{
		quatd test(1, 2, 3, 4);
		REQUIRE(test.scalar() == 1);
		REQUIRE(test.q0() == 1);
		REQUIRE(test.vector() == vec3d(2, 3, 4));
		REQUIRE(test.q1() == 2);
		REQUIRE(test.q2() == 3);
		REQUIRE(test.q3() == 4);

		test = { 0,1,2,3 };
		REQUIRE(test.scalar() == 0);
		REQUIRE(test.q0() == 0);
		REQUIRE(test.vector() == vec3d(1, 2, 3));
		REQUIRE(test.q1() == 1);
		REQUIRE(test.q2() == 2);
		REQUIRE(test.q3() == 3);
	}

	SECTION("Testing fromAxisAndAngle() and getAxisAndAngle()")
	{
		vec3d axis(1, 1, 1);
		axis.normalize();

		const double angle = 0.1;
		q = quatd::fromAxisAndAngle(axis, angle);
		REQUIRE(q.norm() == 1);

		vec3d axis_out;
		double angle_out;
		std::tie(axis_out, angle_out) = q.getAxisAndAngle();
		REQUIRE((relativeError<double>(angle, angle_out) < 2e-14));
		REQUIRE((axis_out - axis).norm() < 2e-14);
	}

	SECTION("Testing slerp()")
	{
		vec3d b(1.2, 1.99, 3.27);
		b.normalize();
		quatd q2 = quatd::fromAxisAndAngle(b, 0.6);

		REQUIRE(quatd::slerp(q, q2, 0.5) == (q*quatd::pow(q.inverse()*q2, 0.5)));
	}

	SECTION("Testing exp() and transform()")
	{
		vec3d x(0, 1, 0);
		const double pi = 3.141592653589793238462643383;
		const vec3d omega(2 * pi, 0, 0);
		const double sec = 10;
		const int N = 100;
		const double dt = sec / N;

		quatd q0, qr;
		for (int i = 0; i < N; i++) {
			quatd integral(0, omega.x()*dt, omega.y()*dt, omega.z()*dt);
			qr = q0*quatd::exp(0.5*integral);
			qr.normalize();

			x = qr.transform(x);
		}

		REQUIRE(x.x() == 0);
		REQUIRE(relativeError<double>(1, x.y()) < 2e-14);
		REQUIRE(x.z() < 2e-15);
	}
}