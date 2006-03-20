/*-
 * Copyright (c) 2006 Ted Stodgell. All rights reserved.
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
 * $Id: Lambert.cpp,v 1.6 2006/03/20 04:10:17 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 *
 */

#include <math.h>
#include "Vector.h"
#include "Global.h"
#include "Lambert.h"
#include <iostream>
using namespace std;

Lambert::Lambert(Vector r1in, Vector r2in, double tin)
{
    double tof = tin;     // time of flight
    Vector r1 = r1in;     // initial vector
    Vector r2 = r2in;     // final vector

    double rr1 = norm(r1);
    double rr2 = norm(r2);
    Vector c12 = cross(r1,r2);
    double theta = acos(dot(r1, r2)/rr1/rr2);
    if ( 0 >= c12.getZ() ) theta = 2 * PI - theta;
    double A = sin(theta) * sqrt(rr1*rr2 / (1-cos(theta)));
}

Lambert::~Lambert (void)
{
    // BE SURE TO FREE DYNAMIC STUFF
}

void
Lambert::psolve (void)
{
    //
}

void
Lambert::rsolve (void)
{
    //
}

double
Lambert::y (double zin)
{
    return 0.0;
}

double
Lambert::F (double zin, double tin)
{
    return 0.0;
}

double
Lambert::dFdz (double zin)
{
    return 0.0;
}


int
main(void) {
    cout << "Testing the Lambert solver.\n";
    Vector p1(1.0, 2.0, 3.0);
    Vector p2(-5.5, -2.1, 0.0);
    double time = 10.0;

    // lambert(p1, p2, time);

    return 0;
}
