/*-
* Copyright (c) 2005 Ted Stodgell. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
* OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
* SUCH DAMAGE.
*
* $Id: Vec3.cpp,v 1.14 2006/08/16 23:36:42 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

#include "Vec3.h"
#include "Orbgnosis.h" // defines DBL_MAX from <float.h>
#include <math.h>
#include <iostream>
using namespace std;

/**
 * The Vec3 constructor with no arguments defaults to all zeroes.
 */
Vec3::Vec3 (void)
{
    for (int i = 0; i < 3; i++)
        e[i] = 0.0;
    // cout << "Default Vec3 constructor called.\n";
}

/**
 * The Vec3 constructor takes three doubles, the vector elements.
 * @param xin is the first element.
 * @param yin is the second element.
 * @param zin is the third element.
 */
Vec3::Vec3 (double xin, double yin, double zin)
{
    e[0] = xin;
    e[1] = yin;
    e[2] = zin;
    // cout << "3-arg Vec3 constructor called.\n";
}

/**
 * The Vec3 destructor.
 */
Vec3::~Vec3 (void)
{
    //cout << "Vec3 dtor.\n";
}

/**
 * The Vec3 copy constructor.
 * @param copy is a reference to the Vec3 to be copied.
 */
Vec3::Vec3 (const Vec3& copy)
{
    for (int i = 0; i < 3; i++)
        e[i] = copy.e[i];
    //cout << "Vec3 copy constructor called\n";
}

/**
 * The Vec3 copy assignment operator.
 */
Vec3&
Vec3::operator = (const Vec3& q)
{
    if (this != &q) // avoid pointless self-assignment
    {
        for (int i = 0; i < 3; i++)
            e[i] = q.e[i];
    }
    return *this;
}

/**
 * The Vec3 addition-assignment operator.
 */
Vec3&
Vec3::operator += (const Vec3& q)
{
    for (int i = 0; i < 3; i++)
        e[i] += q.e[i];
    return *this;
}

/**
 * The Vec3 addition-assignment operator.
 */
Vec3&
Vec3::operator -= (const Vec3& q)
{
    for (int i = 0; i < 3; i++)
        e[i] -= q.e[i];
    return *this;
}

/**
 * Multiplies a Vec3 with a scalar of type double and returns a Vec3.
 * @param q is a vector (type must be Vec3).
 * @param s is a scalar (type must be double).
 */
Vec3
operator * (const Vec3& q, const double& s)
{
    return Vec3 (q.e[0] * s,
                 q.e[1] * s,
                 q.e[2] * s);
}

/**
 * Multiplies a scalar of type double with a Vec3 and returns a Vec3.
 * @param s is a scalar (type must be double).
 * @param q is a vector (type must be Vec3).
 */
Vec3
operator * (const double& s, const Vec3& q)
{
    return Vec3 (q.e[0] * s,
                 q.e[1] * s,
                 q.e[2] * s);
}

/**
 * Divides a Vec3 by a scalar of type double and returns a Vec3.
 * @param q is a vector (type must be Vec3).
 * @param s is a scalar (type must be double).
 */
Vec3
operator / (const Vec3& q, const double& s)
{
    return Vec3 (q.e[0] / s,
                 q.e[1] / s,
                 q.e[2] / s);
}

/**
 * Divides scalar of type double by a Vec3 and returns a Vec3.
 * @param s is a scalar (type must be double).
 * @param q is a vector (type must be Vec3).
 */
Vec3
operator / (const double& s, const Vec3& q)
{
    return Vec3 (s / q.e[0],
                 s / q.e[1],
                 s / q.e[2]);
}

/**
 * Vec3 addition.
 * @param a is a constant reference to the first Vec3.
 * @param b is a constant reference to the second Vec3.
 */
Vec3
operator + (const Vec3& a, const Vec3& b)
{
    return Vec3 (a.e[0] + b.e[0],
                 a.e[1] + b.e[1],
                 a.e[2] + b.e[2]);
}

/**
 * Vec3 unary addition.
 * example: a = + b; 
 */
Vec3
operator + (const Vec3& a)
{
    return a;
}

/**
 * Vec3 subtraction.
 * @param a is a constant reference to the first Vec3.
 * @param b is a constant reference to the second Vec3.
 */
Vec3
operator - (const Vec3& a, const Vec3& b)
{
    return Vec3 (a.e[0] - b.e[0],
                 a.e[1] - b.e[1],
                 a.e[2] - b.e[2]);
}

/**
 * Vec3 unary subtraction.
 */
Vec3
operator - (const Vec3& a)
{
    return Vec3() - a;  // 0 - a = -a
}

/**
 * Cross product, returns type Vec3.
 * @param a is a constant reference to the first Vec3.
 * @param b is a constant reference to the second Vec3.
 */
Vec3
cross (const Vec3& a, const Vec3& b)
{
    const double x = a.e[1] * b.e[2] - a.e[2] * b.e[1];
    const double y = a.e[2] * b.e[0] - a.e[0] * b.e[2];
    const double z = a.e[0] * b.e[1] - a.e[1] * b.e[0];
    return Vec3(x, y, z);
}

/**
 * Dot product, returns type double.
 * @param a is a constant reference to the first Vec3.
 * @param b is a constant reference to the second Vec3.
 */
double
dot (const Vec3& a, const Vec3& b)
{
    return (a.e[0] * b.e[0]) + (a.e[1] * b.e[1]) + (a.e[2] * b.e[2]);
}

/**
 * The Euclidean vector norm.
 * @param q is a constant reference to a Vec3.
 */
double
norm (const Vec3& q)
{
    return sqrt (q.e[0]*q.e[0] + q.e[1]*q.e[1] + q.e[2]*q.e[2]);
}

/**
 * The iostream output operator is overloaded for the Vec3 type.
 * @param s is a reference to ostream.
 * @param q is the Vec3 to be printed with nice formatting. 
 */
ostream&
operator << (ostream& s, Vec3 q)
{
    s << "(" << q.e[0] << ", " << q.e[1] << ", " << q.e[2] << ")";
    return s;
}

/**
 * The iostream input operator is overloaded for the Vec3 type.
 * @param s is a reference to istream.
 * @param q is the Vec3 which we're going to try to read from input.
 */
istream&
operator >> (istream& s, Vec3 q)
{
    s >> q.e[0] >> q.e[1] >> q.e[2];
    return s;
}

/**
 * Rotate a vector about first (x) axis by some angle.
 * It's just a rotation matrix,
 */
Vec3
rotX (const Vec3& v, const double& a)
{
    const double c = cos(a);
    const double s = sin(a);
    // Written out the long way...
    // const double x = 1 * v.e[0] + 0 * v.e[1] + 0 * v.e[2];
    // const double y = 0 * v.e[0] + c * v.e[1] - s * v.e[2];
    // const double z = 0 * v.e[0] + s * v.e[1] + c * v.e[2];
    // or more concisely...
    const double x = v.e[0];
    const double y = c * v.e[1] - s * v.e[2];
    const double z = s * v.e[1] + c * v.e[2];
    return Vec3(x, y, z);
}

/**
 * Rotate a vector about second (y) axis by some angle.
 *  It's just a rotation matrix,
 */
Vec3
rotY (const Vec3& v, const double& a)
{
    const double c = cos(a);
    const double s = sin(a);
    // Written out the long way...
    // const double x =  c * v.e[0] + 0 * v.e[1] + s * v.e[2];
    // const double y =  0 * v.e[0] + 1 * v.e[1] - 0 * v.e[2];
    // const double z = -s * v.e[0] + 0 * v.e[1] + c * v.e[2];
    // or more concisely...
    const double x =  c * v.e[0] + s * v.e[2];
    const double y =  v.e[1];
    const double z = -s * v.e[0] + c * v.e[2];
    return Vec3(x, y, z);
}

/**
 * Rotate a vector about third (z) axis by some angle.
 * It's just a rotation matrix,
 */
Vec3
rotZ (const Vec3& v, const double& a)
{
    const double c = cos(a);
    const double s = sin(a);
    // Written out the long way...
    // const double x = c * v.e[0] - s * v.e[1] + 0 * v.e[2];
    // const double y = s * v.e[0] + c * v.e[1] + 0 * v.e[2];
    // const double z = 0 * v.e[0] + 0 * v.e[1] + 1 * v.e[2];
    // or more concisely...
    const double x = c * v.e[0] - s * v.e[1];
    const double y = s * v.e[0] + c * v.e[1];
    const double z = v.e[2];
    return Vec3(x, y, z);
}

/**
 * Get method for the first element of a Vec3.
 */
double
Vec3::getX (void)
{
    return e[0];
}

/**
 * Get method for the second element of a Vec3.
 */
double
Vec3::getY (void)
{
    return e[1];
}

/**
 * Get method for the third element of a Vec3.
 */
double
Vec3::getZ (void)
{
    return e[2];
}

/**
 * Sets all three elements of a Vec3 to zero.
 */
void
Vec3::toZero (void)
{
    e[0] = e[1] = e[2] = 0.0;
}

/**
 * Sets all three elements of a Vec3 to DBL_MAX (may cause problems).
 */
void
Vec3::toInf (void)
{
    try
    {
        e[0] = e[1] = e[2] = INF;    // INFINITY defined in math.h causes problems
        // with gcc4.0 on Mac OS X.
    }
    catch (...)
    {
        cerr << "Vec3::toInf could not set double to " << INF << "\n";
    }
    return ;
}

/**
 * Sets all three elements of a Vec3 by taking three arguments.
 * @param a is the first element.
 * @param b is the second element.
 * @param c is the third element.
 */
void
Vec3::set3 (double a, double b, double c)
{
    e[0] = a;
    e[1] = b;
    e[2] = c;
}

/**
 * Set method for the first element of a Vec3.
 */
void
Vec3::setX (double xin)
{
    e[0] = xin;
}

/**
 * Set method for the second element of a Vec3.
 */
void
Vec3::setY (double yin)
{
    e[1] = yin;
}

/**
 * Set method for the third element of a Vec3.
 */
void
Vec3::setZ (double zin)
{
    e[2] = zin;
}
