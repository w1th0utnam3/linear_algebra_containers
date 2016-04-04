/*
	Main testing file for the linear_algebra_containers
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

#include <iostream>
#include <string>
#include <cassert>

#include "../source/matrix.h"
#include "../source/column_vector.h"
#include "../source/vector3.h"
#include "../source/quaternion.h"

void msg(const std::string& string) {
	std::cout << string << "... ";
}

void ok() {
	std::cout << "ok." << std::endl;
}

template<typename T>
static T relativeError(T value, T approxValue) {
	using std::abs;
	return abs(T(1)-(value/approxValue));
}

int run_matrix_test()
{
	typedef matrix<double,4,4> mat4x4d;
	typedef matrix<double,4,1> mat4x1d;
	typedef matrix<double,4,2> mat4x2d;

	msg("Testing braced initialization and subscript operator");
	mat4x4d mat{{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}};
	for(int i = 0; i < 16; i++) {
		assert(mat[i] == i+1);
	}
	ok();

	msg("Testing ()-operator");
	assert(mat(0,0) == 1);
	assert(mat(1,0) == 2);
	assert(mat(2,0) == 3);
	assert(mat(3,0) == 4);
	assert(mat(0,1) == 5);
	assert(mat(1,1) == 6);
	assert(mat(2,1) == 7);
	assert(mat(3,1) == 8);
	assert(mat(0,2) == 9);
	assert(mat(1,2) == 10);
	assert(mat(2,2) == 11);
	assert(mat(3,2) == 12);
	assert(mat(0,3) == 13);
	assert(mat(1,3) == 14);
	assert(mat(2,3) == 15);
	assert(mat(3,3) == 16);
	ok();

	msg("Testing fill()");
	mat4x1d v;
	v.fill(12.345);
	assert(v[0] == 12.345);
	assert(v[1] == 12.345);
	assert(v[2] == 12.345);
	assert(v[3] == 12.345);
	ok();

	msg("Testing zeros()");
	v.zeros();
	assert(v[0] == 0);
	assert(v[1] == 0);
	assert(v[2] == 0);
	assert(v[3] == 0);
	ok();

	msg("Testing identity matrix eye()");
	mat.toIdentity();
	for(int i = 0; i < 4; i++) {
		for(int j = 0; j < 4; j++) {
			if(i == j) {
				assert(mat(i,j) == 1);
			} else {
				assert(mat(i,j) == 0);
			}
		}
	}
	ok();

	msg("Testing matrix multiplication");
	for(int i = 0; i < 16; i++) {
		mat[i] = i+1;
	}
	mat4x2d mat2;
	for(int i = 0; i < 8; i++) {
		mat2[i] = i+1;
	}
	matrix<double,4,2> mat3 = mat*mat2;
	assert(mat3(0,0) == 90);
	assert(mat3(1,0) == 100);
	assert(mat3(2,0) == 110);
	assert(mat3(3,0) == 120);
	assert(mat3(0,1) == 202);
	assert(mat3(1,1) == 228);
	assert(mat3(2,1) == 254);
	assert(mat3(3,1) == 280);
	ok();

	msg("Testing matrix scaling");
	matrix<double,2,2> mat4;
	mat4.fill(3);
	mat4 *= 2;
	assert(mat4[0] == 6);
	assert(mat4[1] == 6);
	assert(mat4[2] == 6);
	assert(mat4[3] == 6);
	ok();

	msg("Testing matrix sum/difference");
	matrix<double,2,2> mat5;
	mat5.fill(2);
	mat4 -= mat5;
	assert(mat4[0] == 4);
	assert(mat4[1] == 4);
	assert(mat4[2] == 4);
	assert(mat4[3] == 4);
	mat5.fill(3);
	mat5 += mat4;
	assert(mat5[0] == 7);
	assert(mat5[1] == 7);
	assert(mat5[2] == 7);
	assert(mat5[3] == 7);

	mat4.toIdentity();
	matrix<double,2,2> mat6 = mat5-mat4;
	assert(mat6(0,0) == 6);
	assert(mat6(1,0) == 7);
	assert(mat6(0,1) == 7);
	assert(mat6(1,1) == 6);
	ok();

	msg("Testing matrix comparison");
	assert(mat5 != mat4);
	assert(!(mat5 == mat4));
	ok();

	// TODO: Test data
	// TODO: Test const
	// TODO: Test transposed
	// TODO: Test scalar product

	return 0;
}

int run_column_vector_test()
{
	typedef column_vector<double,3> vec3d;
	vec3d v1;
	vec3d v2;

	// TODO: Test initializer list

	msg("Testing fill() and array subscript operator");
	v1.fill(2.5);
	assert(v1[0] == 2.5);
	assert(v1[1] == 2.5);
	assert(v1[2] == 2.5);
	ok();

	msg("Testing zeors()");
	v2.zeros();
	assert(v2[0] == 0);
	assert(v2[1] == 0);
	assert(v2[2] == 0);
	ok();

	msg("Testing length()");
	v1.fill(3);
	assert(v1.norm() == sqrt(27));
	ok();

	msg("Testing normalize()");
	v1.fill(3);
	v1.normalize();
	assert(v1.norm() == 1);
	ok();

	msg("Testing normalized()");
	v1.fill(3);
	assert(v1.normalized().norm() == 1);
	assert(vec3d::normalized(v1).norm() == 1);
	ok();

	msg("Testing dot product");
	v1[0] = 1;
	v1[1] = 2;
	v1[2] = 3;
	v2[0] = 4;
	v2[1] = 5;
	v2[2] = 6;
	assert(vec3d::dotProduct(v1,v2) == 32);
	ok();

	msg("Testing normSquared");
	v1.fill(2);
	assert(v1.normSquared() == 12);
	v2.fill(3);
	assert(vec3d::normSquared(v2) == 27);
	ok();

	msg("Testing implicit conversion between matrix and column vector");
	v1.fill(2); v2.fill(3);
	assert(vec3d::dotProduct((v1.transposed().transposed()), v2) == 18);
	auto t = v1.transposed().transposed();
	assert(vec3d::dotProduct(t, v2) == 18);
	ok();

	return 0;
}

int run_vector3_test()
{
	typedef vector3<double> vec3d;

	msg("Testing constructor");
	vec3d v1(0.1,312.112,77);
	assert(v1[0] == 0.1);
	assert(v1[1] == 312.112);
	assert(v1[2] == 77);

	const double k = 77;
	vec3d v2(0.1,312.112,k);
	assert(v2[0] == 0.1);
	assert(v2[1] == 312.112);
	assert(v2[2] == k);
	ok();

	msg("Testing getters");
	v1[0] = 3.3;
	v1[1] = 4.4;
	v1[2] = 5.5;
	assert(v1.x() == 3.3);
	assert(v1.y() == 4.4);
	assert(v1.z() == 5.5);
	ok();

	msg("Testing setters");
	const double d = 22.2;
	v2.setX(1);
	v2.setY(d);
	v2.setZ(22);
	assert(v2.x() == 1);
	assert(v2.y() == d);
	assert(v2.z() == 22);
	ok();

	msg("Testing cross product");
	v1.setX(1);
	v1.setY(2);
	v1.setZ(3);
	v2.setX(3);
	v2.setY(4);
	v2.setZ(5);
	auto r = vec3d::crossProduct(v1,v2);
	assert(r.x() == -2);
	assert(r.y() == 4);
	assert(r.z() == -2);
	ok();

	// TODO: Test implicit conversion

	return 0;
}

int run_quaternion_test()
{
	typedef quaternion<double> quatd;
	typedef vector3<double> vec3d;

	msg("Testing constructor");
	vec3d axis(1,1,1);
	axis.normalize();

	const double angle = 0.1;
	quatd q = quatd::fromAxisAndAngle(axis, angle);
	assert(q.norm() == 1);
	ok();

	msg("Testing getAxisAndAngle");
	{

		vec3d axis_out;
		double angle_out;
		std::tie(axis_out, angle_out) = q.getAxisAndAngle();
		assert(relativeError<double>(angle, angle_out) < 2e-14);
		assert(vec3d::norm(axis_out-axis) < 2e-14);
	}
	ok();

	msg("Testing slerp");
	{
		vec3d b(1.2,1.99,3.27);
		b.normalize();
		quatd q2 = quatd::fromAxisAndAngle(b, 0.6);

		assert(quatd::slerp(q,q2,0.5) == (q*quatd::pow(q.inverse()*q2,0.5)));
	}
	ok();

	msg("Testing exp and transform");
	{
		vec3d x(0,1,0);
		const double pi = 3.141592653589793238462643383;
		const vec3d omega(2*pi,0,0);
		const double sec = 10;
		const int N = 100;
		const double dt = sec/N;

		quatd q0, qr;
		for(int i = 0; i < N; i++) {
			quatd integral(0, omega.x()*dt, omega.y()*dt, omega.z()*dt);
			qr = q0*quatd::exp(0.5*integral);
			qr.normalize();

			x = qr.transform(x);
		}

		assert(x.x() == 0);
		assert(relativeError<double>(1, x.y()) < 2e-14);
		assert(x.z() < 2e-15);
	}
	ok();

	// TODO: All tests

	return 0;
}

int main()
{
	#ifdef NDEBUG
		std::cout << "WARNING: This is no debug build! Asserts (i.e. the tests) won't work!" << std::endl;
		std::cout << "Please rebuild without NDEBUG!" << std::endl;
		return 1;
	#endif

	int state = 0;

	std::cout << "Running matrix tests..." << std::endl;
	state += run_matrix_test();
	std::cout << std::endl;

	std::cout << "Running column_vector tests..." << std::endl;
	state += run_column_vector_test();
	std::cout << std::endl;

	std::cout << "Running column_vector tests..." << std::endl;
	state += run_vector3_test();
	std::cout << std::endl;

	std::cout << "Running quaternion tests..." << std::endl;
	state += run_quaternion_test();
	std::cout << std::endl;

	return state;
}
