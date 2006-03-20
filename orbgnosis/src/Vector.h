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
 * $Id: Vector.h,v 1.6 2006/03/20 01:55:22 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 */

#ifndef _VECTOR_H_
#define _VECTOR_H_

class Vector
{
    public:
                    Vector (void);          // no args - 0,0,0

                    Vector (double xin,
                            double yin,
                            double zin);

        virtual     ~Vector (void);

                    Vector (const Vector&); // copy constructor

                    Vector& operator =  (Vector);
                    Vector& operator += (Vector);
                    Vector& operator -= (Vector);

                    double norm (void);
                    void print (void);  // Print vector to stdout.

        friend      Vector cross (const Vector&, const Vector&);
        friend      double dot   (const Vector&, const Vector&);

                    double getX (void);
                    double getY (void);
                    double getZ (void);

    private:
        double x, y, z;
};

Vector operator + (Vector);
Vector operator - (Vector);
Vector operator + (Vector, Vector);
Vector operator - (Vector, Vector);

#endif /* _VECTOR_H_ */
