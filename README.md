# linear_algebra_containers
Some basic, c++ template linear algebra classes (vector, matrix, quaternion)
useful for computer graphics or other computations. Note: the headers requiere
c++11 support to compile.

All files of this project are licensed under the GNU GPL v2 license. Have a look
at the LICENSE file for more information.

## Overview
Classes already implemented:
 - `matrix`: m x n matrix with entries of type T
 - `column_vector`: column vector with n entries of type T
 - `vector3`: column vector with 3 entries of type T
 - `quaternion`: class for rotations etc., with 4 entries of type T

All classes are using templates. For example the `matrix` template parameters
are:
```c++
template<typename T, size_t row_count, size_t column_count>
class matrix;
```
The other classes use an appropriate subset of these parameters as all of them
(except for the quaternion) inherit from the `matrix` class.

## Examples
 Matrix operations are very easy to use and thanks to the template design only
 mathematically correct operations are possible. Example for a matrix product:
 ```c++
 matrix<double,2,3> mat1; mat1.toIdentity();
 matrix<double,3,1> mat2; mat2.toIdentity();
 matrix<double,2,1> result = mat1*mat2;
 ```
The multiplication in reverse order will not work because the matrix dimensions
do not match.
```c++
auto result = mat2*mat1;  // Won't compile
```
Suppose you want to calculate the 2-norm of a column vector. There are three
ways to get the same result:
```c++
typedef column_vector<double,4> vec4;

matrix<double,4,1> x{{1,2,3,4}};
double res1 = std::sqrt(x.transposed()*x);          // 1
double res2 = std::sqrt(vec4::dotProduct(x,x));     // 2
double res3 = vec4(x).length();                     // 3
assert(res1 == res2);
assert(res2 == res3);
```
 - *1* works because a [1 x m] * [m x 1] matrix multiplication returns a
   value of type T instead of a [1x1] matrix.
 - *2* works because a [m x 1] matrix can be implicitly converted to a
   m-dimensional `column_vector`.
 - *3* works because of *2* and thanks to the convenience method of the
   `column_vector`.

All methods are (in my opinion) appropriately documented with doxygen
compatible comments. I will shortly add a doxygen script to generate a
documentation for the classes.

## Todo
This Readme will be updated shortly. There is still much work to do on the
classes even though most basic operations are working. If you want to help or
if you found any bugs feel free to contact me.

Main areas of work:
 - more quaternion methods
 - more tests
 - Matrix4x4 and Matrix3x3 classes with convenience methods
 - more graphics related features like rotations etc.
