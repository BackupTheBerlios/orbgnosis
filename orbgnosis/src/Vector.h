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
 * $Id: Vector.h,v 1.18 2006/06/06 20:45:54 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 */

#ifndef _VECTOR_H_
#define _VECTOR_H_
#include <iostream>
using namespace std;

/**
 * A Vector of exactly three doubles, with operators and common functions.
 * The Vector class overloads vector addition and subtraction,
 * as well as vector-scalar multiplication and division.  Member functions
 * include Euclidean norm, dot product and cross product.  For nicer I/O
 * the >> and << operators are overloaded as well.
 */

class Vector
{
    public:
                    Vector (void);          // defaults to (0,0,0).

                    Vector (double xin,
                            double yin,
                            double zin);

        virtual     ~Vector (void);

                    Vector (const Vector&); // copy constructor

                    Vector& operator =  (Vector);

        // Multiplcation and division with scalar doubles
        // The order of args can be either way.
        friend      Vector  operator * (const Vector&, const double&);
        friend      Vector  operator * (const double&, const Vector&);
        friend      Vector  operator / (const Vector&, const double&);
        friend      Vector  operator / (const double&, const Vector&);

        // Vector æddition and subtraction
        friend      Vector  operator + (const Vector&, const Vector&);
        friend      Vector  operator - (const Vector&, const Vector&);

        // Cross- and dot-products 
        friend      Vector cross (const Vector&, const Vector&);
        friend      double dot   (const Vector&, const Vector&);
        friend      double norm  (const Vector&);

        // Formatted input and output by overloading >> and <<
        friend      ostream& operator << (ostream&, Vector);
        friend      istream& operator >> (istream&, Vector);

        // Standard get-n-set methods.
                    double getX (void);
                    double getY (void);
                    double getZ (void);
                    void toZero (void); // set all elements to zero.
                    void toInf  (void); // set all elements to Infinity.
                    void set3   (double, double, double);
                    void setX   (double);
                    void setY   (double);
                    void setZ   (double);

    private:
        double x; //!< first vector element
        double y; //!< second vector element
        double z; //!< third vector element
};

#endif /* _VECTOR_H_ */
