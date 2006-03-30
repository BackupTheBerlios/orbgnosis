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
 * $Id: Lambert-backup.cpp,v 1.1 2006/03/30 22:52:10 trs137 Exp $
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

    tol = 0.00000001; // error tolerance
    maxloops = 5000;  // iteration limit

    // Magnitudes of vectors
    rr1 = norm(r1);
    rr2 = norm(r2);

    c12 = cross(r1,r2);

    // Swept true anomaly angle between vectors
    theta = acos(dot(r1, r2)/rr1/rr2);

    ratio = 1.0;

    // Initialize the remaining variables to zero for now.
    A = z = f = g = gdot = 0.0;
    v1.setX(0.0);
    v1.setY(0.0);
    v1.setZ(0.0);
    v2.setX(0.0);
    v2.setY(0.0);
    v2.setZ(0.0);

    // theta for prograde solution.
    if ( 0 >= c12.getZ() ) theta = 2 * PI - theta;

    // theta for retrograde solutions.
    // if ( 0 <= c12.getZ() ) theta = 2 * PI - theta;

    A = sin(theta) * sqrt(rr1*rr2 / (1-cos(theta)));

    // Find where F(z, tof) changes sign and use that value
    // for an initial guess.

    z = -100.0;
    while (0 > F(z, tof)) z = z + 0.1;

     cout << "Initial z guess is " << z << "\n";

    // Newton-Raphson iteration
    i = 0;
    while ( (fabs(ratio) > tol) && (i < maxloops))
    {
        i = i + 1;
        ratio = F(z, tof) / dFdz(z);
        z = z - ratio;
    }

    if (i >= maxloops)
        cout << "** WARNING ** max number of iterations exceeded.\n";

    // cout << "After " << i << " iterations z = " << z << "\n";

    f    = 1 - y(z) / rr1;         // Lagrange coeff.
    g    = A * sqrt( y(z) / MU );  // Lagrange coeff.
    gdot = 1 - y(z) / rr2;

    v1 = 1 / g * (r1 - f*r2);
    v2 = 1 / g * (gdot * r2 - r1);

    cout << "Solution: \n\n";

    cout << "v1 (km/s) = " << v1 << "\n";
    cout << "v2 (km/s) = " << v2 << "\n";
}

double
Lambert::y (double zin)
{
    double S = stump_S(zin);
    double C = stump_C(zin);
    double temp = rr1 + rr2 + A * ( zin * S - 1)/sqrt(C);
/*  CHECKS OUT OK
    cout << "-- inside function y \n";
    cout << "-- z           = " << zin << "\n";
    cout << "-- S(z)        = " << S << "\n";
    cout << "-- C(z)        = " << C << "\n";
    cout << "-- A           = " << A << "\n";
    cout << "-- rr1         = " << rr1 << "\n"; 
    cout << "-- rr2         = " << rr2 << "\n";
    cout << "-- y(z)        = " << poop << "\n\n";
*/
    return temp;
}

double
Lambert::F (double zin, double tin)   // For Newton-Raphson iteration
{
    double yy = y(zin);
    double S = stump_S(zin);
    double C = stump_C(zin);
    double yoverc = yy/C;
    
    if (0 >= yy) yy = - yy;

    double temp = pow((yy/C), 1.5) * S
                  + A * sqrt(yy)
                  - ROOTMU * tin;

    cout << "-- inside function F \n";
    cout << "-- z           = " << zin << "\n";
    cout << "-- t           = " << tin << "\n";
    cout << "-- S(z)        = " << S << "\n";
    cout << "-- C(z)        = " << C << "\n";
    cout << "-- A           = " << A << "\n";
    cout << "-- rr1         = " << rr1 << "\n";  
    cout << "-- rr2         = " << rr2 << "\n";
    cout << "-- y(z)        = " << yy << "\n";
    cout << "-- F(z,t)      = " << temp << "\n\n";

    return temp;

    // return pow((yy/C), 1.5) * S 
    //        + A * sqrt(yy) 
    //        - ROOTMU * tin;
}

double
Lambert::dFdz (double zin)           // For Newton-Raphson iteration
{
    if ( 0 == zin)
    {
        double yy = y(0.0);
        return sqrt(2) / 40 * pow(yy, 1.5) 
               + A / 8 * (sqrt(yy) 
               + A * sqrt(1/2/yy));
    }else{
        double yy = y(zin);
        double C = stump_C(zin);
        double S = stump_S(zin);
        return pow((yy/C), 1.5) * (1/2/zin * (C - 3*S/2/C )
               + 3*(S*S)/4/C)
               + A/8*(3*S/C*sqrt(yy)
               + A*sqrt(C/yy));
    }
}

int
main(void) {
    cout << "Testing the Lambert solver.\n";
    Vector p1(5000.0, 10000.0, 2100.0);
    Vector p2(-14600.0, 2500.0, 7000.0);
    double time = 3600.0;

    cout << "Gravitational parameter (km^3/s^2) = " << MU << "\n\n";

    cout << "r1(km) = " << p1 << "\n";
    cout << "r2(km) = " << p2 << "\n";
    cout << "Elapsed time (s) " << time << "\n\n";

    Lambert test(p1, p2, time);
    test.universal();

    return 0;
}
