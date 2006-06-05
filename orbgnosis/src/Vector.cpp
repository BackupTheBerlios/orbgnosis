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
 * $Id: Vector.cpp,v 1.17 2006/06/05 20:51:37 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 */

#include "Vector.h"
#include "Orbgnosis.h" // defines DBL_MAX from <float.h>
#include <math.h>
#include <iostream>
using namespace std;

/*
 * This pragma disables 
 * "remark #981: operands are evaluated in unspecified order"
 * on the Intel C/C++ compiler.  ICC warns about this in unneccessary
 * cases.  Leave NO_ICC_981 undefined to get the warnings.
 */
#ifdef NO_ICC_981
#pragma warning (disable:981)
#endif

/**
 * The Vector constructor with no arguments defaults to all zeroes.
 */
Vector::Vector (void)
{
    x = y = z = 0.0;    // Vectors with no args default to zero.
}

/**
 * The Vector constructor takes three doubles, the vector elements.
 * @param xin is the first element.
 * @param yin is the second element.
 * @param zin is the third element.
 */
Vector::Vector (double xin, double yin, double zin)
{
    x = xin;
    y = yin;
    z = zin;
}

/**
 * The Vector destructor.
 */
Vector::~Vector (void)
{
    //cout << "Vector dtor.\n";
}

/**
 * The Vector copy constructor.
 * @param copy is a reference to the Vector to be copied.
 */
Vector::Vector (const Vector& copy)
{
    x = copy.x;
    y = copy.y;
    z = copy.z;
}

/**
 * The Vector copy assignment operator.
 * @param q is the Vector on the right hand side of the equality.
 */
Vector&
Vector::operator = (Vector q)
{
    x = q.x;
    y = q.y;
    z = q.z;
    return *this;
}

/**
 * Multiplies a Vector with a scalar of type double and returns a Vector.
 * @param q is a vector (type must be Vector).
 * @param s is a scalar (type must be double).
 */
Vector
operator * (const Vector& q, const double& s)
{
    return Vector (q.x * s,
                   q.y * s,
                   q.z * s);
}

/**
 * Multiplies a scalar of type double with a Vector and returns a Vector.
 * @param s is a scalar (type must be double).
 * @param q is a vector (type must be Vector).
 */
Vector
operator * (const double& s, const Vector& q)
{
    return Vector (q.x * s,
                   q.y * s,
                   q.z * s);
}

/**
 * Divides a Vector by a scalar of type double and returns a Vector.
 * @param q is a vector (type must be Vector).
 * @param s is a scalar (type must be double).
 */
Vector
operator / (const Vector& q, const double& s)
{
    return Vector (q.x / s,
                   q.y / s,
                   q.z / s);
}

/**
 * Multiplies scalar of type double by a Vector and returns a Vector.
 * @param s is a scalar (type must be double).
 * @param q is a vector (type must be Vector).
 */
Vector
operator / (const double& s, const Vector& q)
{
    return Vector (q.x / s,
                   q.y / s,
                   q.z / s);
}

/**
 * Vector addition.
 * @param a is a constant reference to the first Vector.
 * @param b is a constant reference to the second Vector.
 */
Vector
operator + (const Vector& a, const Vector& b)
{
    return Vector (a.x + b.x,
                   a.y + b.y,
                   a.z + b.z);
}

/**
 * Vector subtraction.
 * @param a is a constant reference to the first Vector.
 * @param b is a constant reference to the second Vector.
 */
Vector
operator - (const Vector& a, const Vector& b)
{
    return Vector (a.x - b.x,
                   a.y - b.y,
                   a.z - b.z);
}

/**
 * Cross product, returns type Vector.
 * @param a is a constant reference to the first Vector.
 * @param b is a constant reference to the second Vector.
 */
Vector
cross (const Vector& a, const Vector& b)
{
    const double nx = a.y * b.z - b.y * a.z;
    const double ny = a.z * b.x - b.z * a.x;
    const double nz = a.x * b.y - b.x * a.y;
    return Vector(nx, ny, nz);
}

/**
 * Dot product, returns type double.
 * @param a is a constant reference to the first Vector.
 * @param b is a constant reference to the second Vector.
 */
double
dot (const Vector& a, const Vector& b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

/**
 * The Euclidean vector norm.
 * @param q is a constant reference to a Vector.
 */
double
norm (const Vector& q)
{
    return sqrt (q.x*q.x + q.y*q.y + q.z*q.z);
}

/**
 * The iostream output operator is overloaded for the Vector type.
 * @param s is a reference to ostream.
 * @param q is the Vector to be printed with nice formatting. 
 */
ostream&
operator << (ostream& s, Vector q)
{
    s << "(" << q.x << ", " << q.y << ", " << q.z << ")";
    return s;
}

/**
 * The iostream input operator is overloaded for the Vector type.
 * @param s is a reference to istream.
 * @param q is the Vector which we're going to try to read from input.
 */
istream&
operator >> (istream& s, Vector q)
{
    s >> q.x >> q.y >> q.z;
    return s;
}

/**
 * Get method for the first element of a Vector.
 */
double
Vector::getX (void)
{
    return x;
}

/**
 * Get method for the second element of a Vector.
 */
double
Vector::getY (void)
{
    return y;
}

/**
 * Get method for the third element of a Vector.
 */
double
Vector::getZ (void)
{
    return z;
}

/**
 * Sets all three elements of a Vector to zero.
 */
void
Vector::toZero (void)
{
    x = y = z = 0.0;
}

/**
 * Sets all three elements of a Vector to DBL_MAX (may cause problems).
 * NOTE: Using INFINITY as defined in <math.h> works on Linux and FreeBSD,
 * but generates errors with gcc4.0 on Mac OS X because INFINITY (1.0e99)
 * is greater than DBL_MAX.  DBL_MAX is friendlier, but still needs testing.
 * TODO: Add a boolean to the Vector class to track whether a Vector is
 * Very Large.
 */
void
Vector::toInf (void)
{
    try
    {
        x = y = z = DBL_MAX; // INFINITY defined in math.h causes problems
                             // with gcc4.0 on Mac OS X.
    } catch(...) {
        cerr << "Vector::toInf could not set double to <float.h>'s DBL_MAX.\n";
    }
    return;
}

/**
 * Sets all three elements of a Vector by taking three arguments.
 * @param a is the first element.
 * @param b is the second element.
 * @param c is the third element.
 */
void
Vector::set3 (double a, double b, double c)
{
    x = a;
    y = b;
    z = c;
}

/**
 * Set method for the first element of a Vector.
 */
void
Vector::setX (double xin)
{
    x = xin;
}

/**
 * Set method for the second element of a Vector.
 */
void
Vector::setY (double yin)
{
    y = yin;
}

/**
 * Set method for the third element of a Vector.
 */
void
Vector::setZ (double zin)
{
    z = zin;
}

// Re-enable ICC remark #981
#ifdef NO_ICC_981
#pragma warning (default:981)
#endif
