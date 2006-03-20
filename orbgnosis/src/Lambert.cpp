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
 * $Id: Lambert.cpp,v 1.3 2006/03/20 01:55:21 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 *
 * This file contains the first four Stumpff functions C_k(z).
 * Note that z = alpha * chi^2.  
 *
 * Stumpff functions C_k(z) are defined generally by
 * C_k(z) = [1/k!] - [z / (k+2)!] + [z^2 / (k+4)!] + ...
 *
 */

#include <math.h>
#include "Vector.h"
#include <iostream>
using namespace std;

/*
 * Inputs
 *          mu (double, gravitational parameter)
 *          R1 (Vector, start position)
 *          R2 (Vector, final position)
 *          t  (double, time of flight)
 */

void
lambert(double mu, Vector r1, Vector r2, double t)
{
    double rr1 = r1.norm();
    double rr2 = r2.norm();

    cout << "mu is " << mu << ".\n";

    cout << "r1 is ";
    r1.print();
    cout << ".\n";
    
    cout << "r2 is ";
    r2.print();
    cout << ".\n";

    cout << "|r1| is " << rr1 << ".\n";
    cout << "|r2| is " << rr2 << ".\n";

    Vector c12;
    Vector c21;

    cout << "c12 is ";
    c12.print();
    cout << ".\n";

    cout << "c21 is ";
    c21.print();
    cout << ".\n";
}

int
main(void) {
    cout << "Testing the Lambert solver.\n";
    double gravparam = 398600.4418;
    Vector p1(1.0, 2.0, 3.0);
    Vector p2(-5.5, -2.1, 0.0);
    double time = 10.0;

    lambert(gravparam, p1, p2, time);

    return 0;
}
