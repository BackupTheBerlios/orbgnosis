/*
 * Copyright 2006 Ted Stodgell. All rights reserved.
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
 * $Id: Stumpff.cpp,v 1.1 2006/03/19 02:44:10 trs137 Exp $
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
#include <iostream>
using namespace std;

double
stumpff_C0(double z)
{
    if (0.0 == z) return 1.0;        // C0(0) = 1;

    double sqz = 1.0;
    if (z > 0){
        sqz = sqrt(z);
        return ( cos(sqz) );
    }else{ // z < 0
        sqz = sqrt(-z);
        return ( cosh(sqz) );
    }
}

double
stumpff_C1(double z)
{
    if (0.0 == z) return 1.0;        // C1(0) = 1;

    double sqz = 1.0;
    if (z > 0){
        sqz = sqrt(z);
        return ( sin(sqz) / sqz );
    }else{ // z < 0
        sqz = sqrt(-z);
        return ( sinh(sqz) / sqz );
    }
}

double
stumpff_C2(double z)
{
    if (0.0 == z) return 0.5;        // C2(0) = 1/2;

    double sqz = 1.0;
    if (z > 0){
        sqz = sqrt(z);
        return ( (1 - cos(sqz)) / z );
    }else{ // z < 0
        sqz = sqrt(-z);
        return ( (cosh(sqz) - 1) / (-z) );
    }
}

double
stumpff_C3(double z)
{
    if (0.0 == z) return (1.0/6.0);      // C3(0) = 1/6;

    double sqz = 1.0;
    if (z > 0){
        sqz = sqrt(z);
        return ( (sqz - sin(sqz)) / (sqz*sqz*sqz)  );
    }else{ // z < 0
        sqz = sqrt(-z);
        return ( (sinh(sqz) - sqz) /(sqz*sqz*sqz)  );
    }
}

/*
 * To test, uncomment this part and compile Stumpff.cpp alone.
 *
int
main(void) {
    cout << "Stumpff function tester!\n";
    double a = -1.000001;

    while (a < 2.0 )
    {
        cout << "C0(" << a << ") = " << stumpff_C0(a) << "\n";
        cout << "C1(" << a << ") = " << stumpff_C1(a) << "\n";
        cout << "C2(" << a << ") = " << stumpff_C2(a) << "\n";
        cout << "C3(" << a << ") = " << stumpff_C3(a) << "\n\n";
        a = a + 0.20;
    }

    return 0;
}
 */
