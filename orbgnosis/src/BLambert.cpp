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
 * $Id: BLambert.cpp,v 1.2 2006/05/24 15:19:57 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 *                  David Vallado <valladodl@worldnet.att.net>
 *
 */

#include <math.h>
#include "Vector.h"
#include "Global.h"
#include "Stumpff.h"
#include "BLambert.h"
#include <iostream>
using namespace std;

BLambert::BLambert(void)
{
    t = 0.0;
    Ro.toZero();
    R.toZero();
    failure = false;
}

BLambert::BLambert(Vector r1in, Vector r2in, double tin) :
    t(tin), Ro(r1in), R(r2in) 
{
    // cout << "BLambert constructor called \n";
    failure = false;
}

BLambert::~BLambert (void)
{
    // cout << "BLambert destructor called \n";
}

void
BLambert::setRo (Vector vin)
{
    Ro = vin;
}

void
BLambert::setR (Vector vin)
{
    R = vin;
}

void
BLambert::sett (double tin)
{
    t = tin;
}

Vector
BLambert::getVo (void)
{
    return Vo;
}

Vector
BLambert::getV (void)
{
    return V;
}

double
BLambert::gett (void)
{
    return t;
}

bool
BLambert::isFailure(void)
{
    if (true == failure) {
        return true;
    }else{
        return false;
    }
}

/*
 * BATTIN'S METHOD
 *
 * Adapted from David Vallado's implementations
 * "Fundamentals of Astrodynamics and Applications"
 */
void
BLambert::battin (void)
{
    // Local variables
    int     Loops;
    Vector  RCrossR;
    double  u, b, Sinv, Cosv, rp, x, xn, y, L, m, CosDeltaNu,
            SinDeltaNu, DNu, a, tan2w, RoR, h1, h2, tempx, eps,
            denom, chord, k2, s, f, g, fDot, am, ae, be, tm, gDot,
            arg1, arg2, AlpE, AlpH, BetE, BetH, DE, DH, sqdnu;
    double  Ro4, R4;

    // initialize values
    // Magnitudes of Ro and R
    Ro4 = norm(Ro);
    R4  = norm(R);

    CosDeltaNu = dot(Ro, R) / (Ro4 * R4);
    RCrossR    = cross(Ro, R);
    SinDeltaNu = norm(RCrossR) / (Ro4 * R4);
    DNu        = atan2(SinDeltaNu, CosDeltaNu); // quadrant safe

    RoR   = R4 / Ro4;
    eps   = RoR / 1.0;
    tan2w = 0.25 * eps * eps / (sqrt(RoR) + RoR * (2.0 + sqrt (RoR)));
    rp    = sqrt(Ro4 * R4) * ( cos(0.25 * DNu) * cos(0.25 * DNu) + tan2w);

    if (DNu < M_PI)
    {
        sqdnu = sin(0.25 * DNu) * sin (0.25 * DNu);
        L = (sqdnu + tan2w) / (sqdnu + tan2w + cos (0.5 * DNu));
    }else{
        sqdnu = cos(0.25 * DNu) * sin (0.25 * DNu);
        L = (sqdnu + tan2w - cos(0.5 * DNu)) / (sqdnu + tan2w);
    }

    m     = t * t / (8.0 * rp * rp * rp); // t is time of flight
    xn    = 0.0;                          // 0 for parabolic and hyperbolic
    chord = sqrt(Ro4 * Ro4 + R4 * R4 - 2.0 * Ro4 * R4 * cos(DNu));
    s     = 0.5 * (Ro4 + R4 + chord);

    Loops = 1;
    while( Loops < 30)
    {
        x       = xn;
        tempx   = bat_SEE(x);
        denom   = 1.0 / ((1.0 + 2.0 * x + L)
                  * (3.0 + x * (1.0 + 4.0 * tempx)));
        h1      = (L+x)*(L+x) * (1.0 + (1.0 + 3.0 * x) * tempx) * denom;
        h2      = m * (1.0 + (x-L) * tempx) * denom;

        // Evaluate the cubic.
        b   = 0.25 * 27.0 * h2 / ((1.0+h1)*(1.0+h1)*(1.0+h1));
        u   = -0.5 * b / (1.0 + sqrt(1.0 + b));
        k2  = bat_K(u);
        y   = ((1.0 + h1) / 3.0) * (2.0 + sqrt(1.0 + b) /
              (1.0 - 2.0 * u * k2 * k2 ));
        xn  = sqrt( (0.5 * (1.0 - L))*(0.5 * (1.0 - L)) + m / (y*y))
              - 0.5 *  (1.0 + L);

        if (fabs(xn - x) < SMALL) break; // XXX ugly
        Loops = Loops + 1;
    } // end while loop
    a = t * t / (16.0 * rp * rp * xn * y * y);
    // a = rp * m / (2.0 * xn * y * y);  -- XXX commented out in original
    // Find eccentric anomalies
    // Hyperbolic
    if (a < -SMALL)
    {
        arg1 = sqrt(s / (-2.0 * a));
        arg2 = sqrt((s - chord) / (-2.0 * a));

        // evaluate f and g functions

        /*
         * C++ note: asinh() acosh() and atanh() are not in the C90
         * standard.
         * doing #include <math.h> may or may not transparently
         * include the appropriate header, ymmv.
         *
         * asinh() and friends will work, if...
         *
         * FreeBSD: in /usr/include/math.h, only if 
         * #if __BSD_VISIBLE || __ISO_C_VISIBLE >= 1999 || __XSI_VISIBLE
         *
         * Linux RHEL: [TODO] it's not in /usr/include/math.h... yet it works?
         *
         * AIX 4.3: in /usr/include/math.h, only if
         * #if _XOPEN_SOURCE_EXTENDED==1
         *
         * MacOS X: [TODO]
         */
        AlpH = 2.0 * asinh(arg1);

    } // end if

} // end BLambert::battin

double
BLambert::bat_SEE(double v)
{
    // c: array (0..20) of Real;
    double  c[20];
    double  term, termold, del, delold, sum1, temp, eta, SQRTopv;
    int     i;

    c[0] = 0.2;
    for (int j = 1; j<21; j++)
    {
        temp = (double)(j + 2);
        c[j] = (temp * temp) / ((4 * temp * temp) - 1);
    }
    SQRTopv = sqrt(1.0 + v);
    eta = v / (1.0 + 2*SQRTopv + SQRTopv*SQRTopv);

    // Process Forwards
    delold  = 1.0;
    termold = c[0]; // * eta
    sum1    = termold;
    i       = 1;
    while ((i<=20) && (fabs(termold) > 0.000001))
    {
        del     = 1.0 / (1.0 + c[i] * eta * delold);
        term    = termold * (del - 1.0);
        sum1    = sum1 + term;
        i       = i + 1;
        delold  = del;
        termold = term;
    } // end while loop

    return (1.0 / (8.0 * (1.0 + SQRTopv))) *
           (3.0 + sum1 / (1.0 + eta * sum1));
} // end BLambert::bat_SEE

double
BLambert::bat_K(double v)
{
    // d: array (0..20) of Real; -- hardcoded, see astiod.adb
    /*
    double d[21];
    d[0]    = 0; // Padding. Not used.
    d[1]    = 1.0       / 3.0;
    d[2]    = 4.0       / 27.0;
    d[3]    = 8.0       / 27.0;
    d[4]    = 2.0       / 9.0;
    d[5]    = 22.0      / 81.0;
    d[6]    = 208.0     / 891.0;
    d[7]    = 340.0     / 1287.0;
    d[8]    = 418.0     / 1755.0;
    d[9]    = 598.0     / 2295.0;
    d[10]   = 700.0     / 2907.0;
    d[11]   = 928.0     / 3591.0;
    d[12]   = 1054.0    / 4347.0;
    d[13]   = 1330.0    / 5175.0;
    d[14]   = 1480.0    / 6075.0;
    d[15]   = 1804.0    / 7047.0;
    d[16]   = 1978.0    / 8091.0;
    d[17]   = 2350.0    / 9207.0;
    d[18]   = 2548.0    / 10395.0;
    d[19]   = 2968.0    / 11655.0;
    d[20]   = 3190.0    / 12987.0;
    d[21]   = 3658.0    / 14391.0;
    */

    // Static function variables are initialized once and only one
    // copy is created even if the function is called recursively.

    const static double d[] =
    {   0,
        (1.0       / 3.0 ),
        (4.0       / 27.0),
        (8.0       / 27.0),
        (2.0       / 9.0 ),
        (22.0      / 81.0),
        (208.0     / 891.0),
        (340.0     / 1287.0),
        (418.0     / 1755.0),
        (598.0     / 2295.0),
        (700.0     / 2907.0),
        (928.0     / 3591.0),
        (1054.0    / 4347.0),
        (1330.0    / 5175.0),
        (1480.0    / 6075.0),
        (1804.0    / 7047.0),
        (1978.0    / 8091.0),
        (2350.0    / 9207.0),
        (2548.0    / 10395.0),
        (2968.0    / 11655.0),
        (3190.0    / 12987.0),
        (3658.0    / 14391.0)
    };

    int     i;
    double  del, delold, term, termold, sum1;

    // process fowards
    sum1 = d[0]; // XXX check wtf i was thinking
    delold = 1.0;
    termold = d[0];
    i = 1;
    while ((i <= 20) && ( fabs(termold > 0.000001) ))
    {
        del     = 1.0 / (1.0 - d[i] * v * delold);
        term    = termold * (del - 1.0);
        sum1    = sum1 + term;
        i       = i + 1;
        delold  = del;
        termold = term;
    } // end while loop
    return sum1;
} // end BLambert::bat_K
