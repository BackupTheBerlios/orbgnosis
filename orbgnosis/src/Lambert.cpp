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
 * $Id: Lambert.cpp,v 1.13 2006/03/27 18:23:17 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 *
 */

#include <math.h>
#include "Vector.h"
#include "Global.h"
#include "Stumpff.h"
#include "Lambert.h"
#include <iostream>
using namespace std;

/* Generic Lambert constructor to be used for all solution methods */
Lambert::Lambert(const Vector r1in, const Vector r2in, const double tin) :
    tof(tin), r1(r1in), r2(r2in) // Constants. Member initialization list.
{
    cout << "Lambert constructor called \n";
}

Lambert::~Lambert (void)
{
    // BE SURE TO FREE DYNAMIC STUFF
    cout << "Lambert destructor called \n";
}


/*
 * UNIVERSAL VARIABLES METHOD
 *
 * Member functions are:
 *          universal, y, F, dFdz
 */

void
Lambert::universal (void)  // Prograde universal solution
{
    cout << "Beginning universal variable solution.\n";
    // Magnitudes of vectors
    rr1 = norm(r1);
    rr2 = norm(r2);

    c12 = cross(r1,r2);

    // Swept angle
    theta = acos(dot(r1, r2)/rr1/rr2);

    // Initialize the remaining variables to zero for now.
    A = z = f = g = gdot = 0.0;
    v1.setX(0.0);
    v1.setY(0.0);
    v1.setZ(0.0);
    v2.setX(0.0);
    v2.setY(0.0);
    v2.setZ(0.0);

    // PROGRADE
    if ( 0 >= c12.getZ() ) theta = 2 * PI - theta;
    double A = sin(theta) * sqrt(rr1*rr2 / (1-cos(theta)));
}

double
Lambert::y (double zin)
{
    return rr1 + rr2 + A * ( zin * stumpff_C3(zin) - 1)/sqrt(stumpff_C2(zin));
}

double
Lambert::F (double zin, double tin)   // For Newton-Raphson iteration
{
    double my_y = y(zin);
    return pow( my_y/stumpff_C2(zin), 1.5) * stumpff_C3(zin) 
           + A * sqrt(my_y) 
           - ROOTMU * tin;
}

double
Lambert::dFdz (double zin)           // For Newton-Raphson iteration
{
    if ( 0 == zin)
    {
        double my_y = y(0.0);
        return sqrt(2) / 40* pow(my_y, 1.5) 
               + A / 8 * (sqrt(my_y) 
               + A * sqrt(1/2/my_y));
    }else{
        double my_y = y(zin);
        double C = stumpff_C2(zin);
        double S = stumpff_C3(zin);

        return pow(y(z)/C,1.5) * (1/2/zin*(C - 3*S/2/C )
               + 3*S*S/4/C)
               + A/8*(3*S/C*sqrt(my_y)
               + A*sqrt(C/my_y));
    }
}


int
main(void) {
    cout << "Testing the Lambert solver.\n";
    Vector p1(1.0, 2.0, 3.0);
    Vector p2(-5.5, -2.1, 0.0);
    double time = 10.0;

    Lambert test(p1, p2, time);
    test.universal();

    return 0;
}
