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
* $Id: Vec3.cpp,v 1.1 2006/06/15 20:50:33 trs137 Exp $
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
    : x(0.0),
      y(0.0),
      z(0.0)
{
    // cout << "Default Vec3 constructor called.\n";
}

/**
 * The Vec3 constructor takes three doubles, the vector elements.
 * @param xin is the first element.
 * @param yin is the second element.
 * @param zin is the third element.
 */
Vec3::Vec3 (double xin, double yin, double zin)
        : x(xin),
        y(yin),
        z(zin)
{
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
        : x(copy.x),
        y(copy.y),
        z(copy.z)
{
    //cout << "Vec3 copy constructor called\n";
}

/**
 * The Vec3 copy assignment operator.
 * @param q is the Vec3 on the right hand side of the equality.
 */
Vec3&
Vec3::operator = (Vec3 q)
{
    x = q.x;
    y = q.y;
    z = q.z;
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
    return Vec3 (q.x * s,
                   q.y * s,
                   q.z * s);
}

/**
 * Multiplies a scalar of type double with a Vec3 and returns a Vec3.
 * @param s is a scalar (type must be double).
 * @param q is a vector (type must be Vec3).
 */
Vec3
operator * (const double& s, const Vec3& q)
{
    return Vec3 (q.x * s,
                   q.y * s,
                   q.z * s);
}

/**
 * Divides a Vec3 by a scalar of type double and returns a Vec3.
 * @param q is a vector (type must be Vec3).
 * @param s is a scalar (type must be double).
 */
Vec3
operator / (const Vec3& q, const double& s)
{
    return Vec3 (q.x / s,
                   q.y / s,
                   q.z / s);
}

/**
 * Divides scalar of type double by a Vec3 and returns a Vec3.
 * @param s is a scalar (type must be double).
 * @param q is a vector (type must be Vec3).
 */
Vec3
operator / (const double& s, const Vec3& q)
{
    return Vec3 (s / q.x,
                   s / q.y,
                   s / q.z);
}

/**
 * Vec3 addition.
 * @param a is a constant reference to the first Vec3.
 * @param b is a constant reference to the second Vec3.
 */
Vec3
operator + (const Vec3& a, const Vec3& b)
{
    return Vec3 (a.x + b.x,
                   a.y + b.y,
                   a.z + b.z);
}

/**
 * Vec3 subtraction.
 * @param a is a constant reference to the first Vec3.
 * @param b is a constant reference to the second Vec3.
 */
Vec3
operator - (const Vec3& a, const Vec3& b)
{
    return Vec3 (a.x - b.x,
                   a.y - b.y,
                   a.z - b.z);
}

/**
 * Cross product, returns type Vec3.
 * @param a is a constant reference to the first Vec3.
 * @param b is a constant reference to the second Vec3.
 */
Vec3
cross (const Vec3& a, const Vec3& b)
{
    const double nx = a.y * b.z - b.y * a.z;
    const double ny = a.z * b.x - b.z * a.x;
    const double nz = a.x * b.y - b.x * a.y;
    return Vec3(nx, ny, nz);
}

/**
 * Dot product, returns type double.
 * @param a is a constant reference to the first Vec3.
 * @param b is a constant reference to the second Vec3.
 */
double
dot (const Vec3& a, const Vec3& b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

/**
 * The Euclidean vector norm.
 * @param q is a constant reference to a Vec3.
 */
double
norm (const Vec3& q)
{
    return sqrt (q.x*q.x + q.y*q.y + q.z*q.z);
}

/**
 * The iostream output operator is overloaded for the Vec3 type.
 * @param s is a reference to ostream.
 * @param q is the Vec3 to be printed with nice formatting. 
 */
ostream&
operator << (ostream& s, Vec3 q)
{
    s << "(" << q.x << ", " << q.y << ", " << q.z << ")";
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
    s >> q.x >> q.y >> q.z;
    return s;
}

/**
 * Get method for the first element of a Vec3.
 */
double
Vec3::getX (void)
{
    return x;
}

/**
 * Get method for the second element of a Vec3.
 */
double
Vec3::getY (void)
{
    return y;
}

/**
 * Get method for the third element of a Vec3.
 */
double
Vec3::getZ (void)
{
    return z;
}

/**
 * Sets all three elements of a Vec3 to zero.
 */
void
Vec3::toZero (void)
{
    x = y = z = 0.0;
}

/**
 * Sets all three elements of a Vec3 to DBL_MAX (may cause problems).
 * TODO: Add a boolean to the Vec3 class to track whether a Vec3 is
 * Very Large.
 */
void
Vec3::toInf (void)
{
    try
    {
        x = y = z = INF;    // INFINITY defined in math.h causes problems
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
    x = a;
    y = b;
    z = c;
}

/**
 * Set method for the first element of a Vec3.
 */
void
Vec3::setX (double xin)
{
    x = xin;
}

/**
 * Set method for the second element of a Vec3.
 */
void
Vec3::setY (double yin)
{
    y = yin;
}

/**
 * Set method for the third element of a Vec3.
 */
void
Vec3::setZ (double zin)
{
    z = zin;
}
