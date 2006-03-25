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
 * $Id: Vector.cpp,v 1.10 2006/03/25 21:52:34 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 */

#include "Vector.h"
#include <math.h>
#include <iostream>
using namespace std;

Vector::Vector (void)
{
    x = y = z = 0.0;
    cout << "Vector ctor.\n";
}

Vector::Vector (double xin, double yin, double zin)
{
    x = xin;
    y = yin;
    z = zin;
    cout << "Vector ctor.\n";
}

Vector::~Vector (void)
{
    //
    cout << "Vector dtor.\n";
}

Vector::Vector (const Vector& copy)
{
    x = copy.x;
    y = copy.y;
    z = copy.z;
}

Vector&
Vector::operator = (Vector q)
{
    x = q.x;
    y = q.y;
    z = q.z;
    return *this;
}

Vector&
Vector::operator += (Vector q)
{
    x += q.x;
    y += q.y;
    z += q.z;
    return *this;
}

Vector&
Vector::operator -= (Vector q)
{
    x -= q.x;
    y -= q.y;
    z -= q.z;
    return *this;
}

Vector
cross (const Vector &a, const Vector &b)
{
    const double nx = a.y * b.z - b.y * a.z;
    const double ny = a.z * b.x - b.z * a.x;
    const double nz = a.x * b.y - b.x * a.y;
    return Vector(nx, ny, nz);
}

double
dot (const Vector &a, const Vector &b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

double
norm (const Vector&q)
{
    return sqrt (q.x*q.x + q.y*q.y + q.z*q.z);
}

ostream&
operator << (ostream& s, Vector q)
{
    s << "(" << q.x << ", " << q.y << ", " << q.z << ")";
    return s;
}

istream&
operator >> (istream& s, Vector q)
{
    s >> q.x >> q.y >> q.z;
    return s;
}

double
Vector::getX (void)
{
    return x;
}

double
Vector::getY (void)
{
    return y;
}

double
Vector::getZ (void)
{
    return z;
}

void
Vector::setX (double xin)
{
    x = xin;
}

void
Vector::setY (double yin)
{
    y = yin;
}
void
Vector::setZ (double zin)
{
    z = zin;
}

inline Vector
operator + (Vector q)
{
    return q;
}

inline Vector
operator - (Vector q)
{
    return -q;
}

inline Vector
operator + (Vector a, Vector b)
{
    return a += b;
}

inline Vector
operator - (Vector a, Vector b)
{
    return a -= b;
}
