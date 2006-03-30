/*-
 * Copyright 2006 (c) Ted Stodgell. All rights reserved.
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
 * $Id: Lambert.h,v 1.7 2006/03/30 05:36:19 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 */

#ifndef _LAMBERT_H_
#define _LAMBERT_H_

class Lambert 
{
    public:
                Lambert (const Vector r1in,
                         const Vector r2in,
                         const double tin);

                virtual ~Lambert (void);

                Lambert (const Lambert&); // copy ctor

                /* Universal Variable method */
                void universal(void);
                double y (double zin);
                double F(double zin, double tin);
                double dFdz(double zin);
                /*****************************/

                /* Battin's method */
                // TODO
                /*****************************/

    private:
        const double tof;
        const Vector r1;
        const Vector r2;

        double tol;             // error tolerance
        int maxloops;           // iteration limit
        int i;                  // loop index

        double rr1;
        double rr2;
        Vector c12;
        double theta;           // angle between r1 and r2;
        double A;
        double z;
        double f;               // Lagrange coefficient
        double g;               // Lagrange coefficient
        double gdot;            // dg/dt
        double ratio;           // F / dFdz
        Vector v1;              // initial velocity of xfer arc
        Vector v2;              // final velocity of xfer arc
};

#endif /* _LAMBERT_H_ */
