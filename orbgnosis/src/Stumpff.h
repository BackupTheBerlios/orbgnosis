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
 * $Id: Stumpff.h,v 1.7 2006/06/05 20:51:37 trs137 Exp $
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

/** @file
  * Stumpff functions are used in the universal variables solution to
  * the orbital equation.  The Stumpff function is defined by,
  * C_k(z) = 1/k! - x/(k+2)! + x^2 / (k+4)! + ...
  * This header file defines the first four functions, C_0 through C_3.
  */
#ifndef _STUMPFF_H_
#define _STUMPFF_H_

#include <math.h>
#include <iostream>

/**
 * The zeroth Stumpff c function
 */
inline double
stumpff_C0(double z)
{
    if (fabs(z) < SMALL) return 1.0;        // C0(0) = 1;

    double sqz = 1.0;
    if (z > 0){
        sqz = sqrt(z);
        return ( cos(sqz) );
    }else{ // z < 0
        sqz = sqrt(-z);
        return ( cosh(sqz) );
    }
}

/**
 * The first Stumpff c function
 */
inline double
stumpff_C1(double z)
{
    if (fabs(z) < SMALL) return 1.0;        // C1(0) = 1;

    double sqz = 1.0;
    if (z > 0){
        sqz = sqrt(z);
        return ( sin(sqz) / sqz );
    }else{ // z < 0
        sqz = sqrt(-z);
        return ( sinh(sqz) / sqz );
    }
}

/**
 * The second Stumpff c function
 */
inline double
stumpff_C2(double z)
{
    if (fabs(z) < SMALL) return 0.5;        // C2(0) = 1/2;

    double sqz = 1.0;
    if (z > 0){
        sqz = sqrt(z);
        return ( (1 - cos(sqz)) / z );
    }else{ // z < 0
        sqz = sqrt(-z);
        return ( (cosh(sqz) - 1) / (-z) );
    }
}

/**
 * The third Stumpff c function
 */
inline double
stumpff_C3(double z)
{
    if (fabs(z) < SMALL) return (1.0/6.0);      // C3(0) = 1/6;

    double sqz = 1.0;
    if (z > 0){
        sqz = sqrt(z);
        return ( (sqz - sin(sqz)) / (sqz*sqz*sqz)  );
    }else{ // z < 0
        sqz = sqrt(-z);
        return ( (sinh(sqz) - sqz) /(sqz*sqz*sqz)  );
    }
}

#endif // _STUMPFF_H_
