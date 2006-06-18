/*-
* Copyright 2005 (c) Ted Stodgell. All rights reserved.
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
* $Id: Vec3.h,v 1.3 2006/06/18 23:47:18 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

#ifndef _VEC3_H_
#define _VEC3_H_
#include <iostream>
using namespace std;

/**
 * A Vec3 of exactly three doubles, with operators and common functions.
 * The Vec3 class overloads vector addition and subtraction,
 * as well as vector-scalar multiplication and division.  Member functions
 * include Euclidean norm, dot product and cross product.  For nicer I/O
 * the >> and << operators are overloaded as well.
 */
class Vec3
{
    public:
        Vec3 (void);          // defaults to (0,0,0).
        Vec3 (double, double, double);

        virtual ~Vec3 (void);

        Vec3 (const Vec3&);   // copy constructor

        Vec3& operator = (const Vec3&);

        Vec3& operator += (const Vec3&); // add-assign
        Vec3& operator -= (const Vec3&); // subtract-assign

        // Scalar Multiplcation and division.
        // The order of args can be either way.
        friend Vec3 operator * (const Vec3&, const double&);
        friend Vec3 operator * (const double&, const Vec3&);
        friend Vec3 operator / (const Vec3&, const double&);
        friend Vec3 operator / (const double&, const Vec3&);

        // Vec3 addition and subtraction.
        friend Vec3 operator + (const Vec3&, const Vec3&); // binary
        friend Vec3 operator + (const Vec3&);              // unary
        friend Vec3 operator - (const Vec3&, const Vec3&); // binary
        friend Vec3 operator - (const Vec3&);              // unary

        // Cross- and dot-products.
        friend Vec3 cross (const Vec3&, const Vec3&);
        friend double dot (const Vec3&, const Vec3&);
        friend double norm (const Vec3&);

        // Formatted input and output by overloading >> and <<
        friend ostream& operator << (ostream&, Vec3);
        friend istream& operator >> (istream&, Vec3);

        // Standard get-n-set methods.
        double getX (void);
        double getY (void);
        double getZ (void);
        void toZero (void); // set all elements to zero.
        void toInf (void); // set all elements to Infinity.
        void set3 (double, double, double);
        void setX (double);
        void setY (double);
        void setZ (double);

    private:
        double e[3];
};

#endif /* _VEC3_H_ */
