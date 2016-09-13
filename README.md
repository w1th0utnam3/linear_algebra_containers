# linear_algebra_containers
Some basic, c++ template linear algebra classes (vector, matrix, quaternion)
useful for non high performance calculations (the classes don't use template
expressions, unrolling or other optimizations).

All files of this project are licensed under the MIT license. Have a look
at the LICENSE file for more information.

## Overview
Classes already implemented:
 - `Matrix`: m x n matrix with entries of type T
 - `ColumnVector`: column vector with n entries of type T
 - `Vector3`: column vector with 3 entries of type T
 - `Quaternion`: class for rotations etc., with 4 entries of type T

All classes are using templates. For example the `Matrix` template parameters
are:
```c++
template<typename T, size_t row_count, size_t column_count>
class Matrix;
```
The `ColumnVector` is a partial specialization of the corresponding matrix type
while the Vector3 class is a subclass of the ColumnVector.

## Examples
 Matrix operations are very easy to use and thanks to the template design only
 mathematically correct operations are possible. Example for a matrix product:
 ```c++
 lin_algebra::Matrix<double,2,3> mat1; mat1.toIdentity();
 lin_algebra::Matrix<double,3,1> mat2; mat2.toIdentity();
 lin_algebra::Matrix<double,2,1> result = mat1*mat2;
 ```
The multiplication in reverse order will not work because the matrix dimensions
do not match.
```c++
auto result = mat2*mat1;    // Won't compile
```
Suppose you want to calculate the 2-norm of a column vector. There are three
ways to get the same result:
```c++
typedef lin_algebra::ColumnVector<double,4> vec4;
typedef lin_algebra::Matrix<double,4,1> mat4x1;

vec4 x{1.,2.,3.,4.};
mat4x1 y{1.,2.,3.,4.};

double res1 = std::sqrt(x.transposed()*x);          // 1
double res2 = std::sqrt(vec4::dotProduct(x,x));     // 2
double res3 = y.norm();                             // 3
assert(res1 == res2);
assert(res2 == res3);
```
`ColumnVector<T,n>` is just an alias for `Matrix<T,n,1>` so it's not surprising
that everything above works. Additionally the partial specialization allows the
following to work:
```c++
typedef lin_algebra::Matrix<double,4,4> mat4x4;

mat4x4 eye;
eye.toIdentity();
vec4 empty{0;0;0;0};

double res4 = (M*x + empty).norm();
assert(res4 == res1);
```

## Todo
There is still much work to do on the classes even though most basic operations
are working. If you want to help or if you find any bugs feel free to contact
me.

Main areas of work:
 - More quaternion methods
 - More tests
 - Matrix4x4 and Matrix3x3 classes with convenience methods
 - More graphics related features like rotations etc.
 - Point classes & coordinate system transformations
