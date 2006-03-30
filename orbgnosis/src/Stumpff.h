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
 * $Id: Stumpff.h,v 1.2 2006/03/30 05:55:11 trs137 Exp $
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
#ifndef _STUMPFF_H_
#define _STUMPFF_H_

#include <math.h>
#include <iostream>

/*inline*/ double
stump_C(double z)
{
    if (0.0 == z) return 0.5;        // C2(0) = 1/2;

    if (z > 0){
        return ( (1 - cos(sqrt(z))) / z );
    }else{ // z < 0
        return ( (cosh(sqrt(-z)) - 1) / (-z) );
    }
}

/*inline*/ double
stump_S(double z)
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

#endif // _STUMPFF_H_
