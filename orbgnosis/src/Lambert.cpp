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
 * $Id: Lambert.cpp,v 1.37 2006/05/24 14:19:19 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 *                  David Vallado <valladodl@worldnet.att.net>
 *
 */

#include <math.h>
#include "Vector.h"
#include "Global.h"
#include "Stumpff.h"
#include "Lambert.h"
#include <iostream>
using namespace std;

Lambert::Lambert(void)
{
    t = 0.0;
    Ro.toZero();
    R.toZero();
    failure = false;
}

Lambert::Lambert(Vector r1in, Vector r2in, double tin) :
    t(tin), Ro(r1in), R(r2in) 
{
    // cout << "Lambert constructor called \n";
    failure = false;
}

Lambert::~Lambert (void)
{
    // cout << "Lambert destructor called \n";
}

void
Lambert::setRo (Vector vin)
{
    Ro = vin;
}

void
Lambert::setR (Vector vin)
{
    R = vin;
}

void
Lambert::sett (double tin)
{
    t = tin;
}

Vector
Lambert::getVo (void)
{
    return Vo;
}

Vector
Lambert::getV (void)
{
    return V;
}

double
Lambert::gett (void)
{
    return t;
}

bool
Lambert::isFailure(void)
{
    if (true == failure) {
        return true;
    }else{
        return false;
    }
}

/*
 * UNIVERSAL VARIABLES METHOD
 *
 * Adapted from David Vallado's implementations
 * "Fundamentals of Astrodynamics and Applications"
 */
void
Lambert::universal (const bool Lin, const int multirev)
{ 
    // Local variables
    const bool longway = Lin;
    const int NumIter = 40;
    int Loops, YNegKtr;
    double VarA, Y, Upper, Lower, CosDeltaNu, F, G, GDot, XOld,
           XOldCubed, PsiOld, PsiNew, C2New, C3New, dtNew;
    double Ro4, R4;

    const int revs = multirev;

    failure = false;

    PsiNew = 0.0;
    Vo.toZero();
    V.toZero();

    // Magnitudes of Ro and R
    Ro4 = norm(Ro);
    R4  = norm(R);

    // "Nu" is true anomaly.
    CosDeltaNu = dot(Ro, R) / (Ro4 * R4);

    if (true == longway)
    {
        VarA = -sqrt( Ro4 * R4 * (1.0 + CosDeltaNu));
    }else{
        VarA = sqrt( Ro4 * R4 * (1.0 + CosDeltaNu));
    }

    // Form initial guesses.
    PsiOld = 0.0;
    PsiNew = 0.0;
    XOld   = 0.0;
    C2New  = 0.5;
    C3New  = 1.0/6.0;

    // Set up initial bounds for the bisection.
    if (0 == revs)
    {
        Upper = 4.0 * M_PI * M_PI;
        Lower = -8.0 * M_PI;
    }else{
        Upper = -SMALL + 4.0 * ((revs/2)+1)*((revs/2)+1)*M_PI*M_PI;
        Lower = SMALL + 4.0 * (revs/2)*(revs/2) * M_PI*M_PI;
    }

    // Determine if the orbit is possible at all
    if (fabs(VarA) > SMALL)
    {
        Loops = 0;
        YNegKtr = 1;  // y neg counter
        dtNew = -10.0;
        while ((fabs(dtNew - t) > SMALL) && (NumIter > Loops))
        {
            if (fabs(C2New) > SMALL )
            {
                Y = Ro4 + R4 -
                    (VarA * (1.0 - PsiOld * C3New) / sqrt(C2New));
            }else{
                Y = Ro4 + R4;
            }
            // Check for negative values of Y.
            if ( (0 < VarA) && (0 > Y))
            {
                YNegKtr = 1;
                while ((0 > Y)&&(10 > YNegKtr))
                {
                    PsiNew = 0.8 * (1/C3New) * (1.0
                             - (Ro4 + R4) * sqrt(C2New)/VarA);
                    // Find C2 and C3 functions.
                    C2New = stumpff_C2(PsiNew);
                    C3New = stumpff_C3(PsiNew);
                    PsiOld = PsiNew;
                    Lower = PsiOld;
                    if (fabs(C2New) > SMALL)
                    {
                        Y = Ro4 + R4 - (VarA * (1.0 -
                             PsiOld*C3New)/sqrt(C2New) );
                    }else{
                        Y = Ro4 + R4;
                    }
                    YNegKtr = YNegKtr + 1;
                } // end while loop.
            } // end if Y neg.

            if (10 > YNegKtr)
            {
                if (fabs(C2New) > SMALL)
                {
                    XOld = sqrt(Y/C2New);
                }else{
                    XOld = 0.0;
                }
                XOldCubed = XOld * XOld * XOld;
                dtNew = XOldCubed * C3New + VarA*sqrt(Y);

                // Readjust upper and lower bounds.
                if (dtNew < t)
                {
                    Lower = PsiOld;
                }else{
                    Upper = PsiOld;
                }
                PsiNew = (Upper+Lower) * 0.5;

                // Find C2 and C3 functions.
                C2New = stumpff_C2(PsiNew);
                C3New = stumpff_C3(PsiNew);
                PsiOld = PsiNew;
                Loops = Loops + 1;

                /********************************************
                DEBUGGING OUTPUT
                cout << "\nIteration   : " << Loops << "\n";
                cout << "Y(ER)       : " << Y << "\n";
                cout << "Xo(sqrt(ER)): " << XOld << "\n";
                cout << "dtNew(TU)   : " << dtNew << "\n";
                cout << "PsiOld      : " << PsiOld << "\n";
                cout << "PsiNew      : " << PsiNew << "\n";
                *********************************************/

                // Make sure the first guess isn't too close.
                if (( fabs(dtNew - t) < SMALL) && (1 == Loops))
                {
                    dtNew = t - 1.0;
                }
            } // end if (10 > YNegKtr)
        } // end while loop

        if ( (Loops >= NumIter) || (YNegKtr > 10) )
        {
            // cout << "Error: Lambert Universal failed to converge. \n";
            //if (YNegKtr > 10) cout << "Y is negative\n";
            //cout << "NumIter = " << NumIter << "\n";
            Vo.set3(INFINITY, INFINITY, INFINITY);
            V.set3(INFINITY, INFINITY, INFINITY);
            failure = true;
        }else{

            // Use F and G series to find velocity vectors.
            F    = 1.0 - Y/Ro4;
            GDot = 1.0 - Y/R4;
            G    = 1.0 / (VarA * sqrt(Y)); // 1 over G

            Vo = (R - F*Ro) * G;
            V  = (GDot * R - Ro) * G;
        } // end if answer has converged
    }else{
        // cout << "Vectors are 180 degrees apart.\n";
        Vo.set3(INFINITY, INFINITY, INFINITY);
        V.set3(INFINITY, INFINITY, INFINITY);
        failure = true;
    } // end if VarA > SMALL

/* Convert from canonical units back to S.I.
    Vo = Vo * ER / TU_SEC;
    V = V * ER / TU_SEC;
*/
} // end Lambert::universal

/*
 * BATTIN'S METHOD
 *
 * Adapted from David Vallado's implementations
 * "Fundamentals of Astrodynamics and Applications"
 */
void
Lambert::battin (void)
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

        if (fabs(xn -x) < SMALL) break; // XXX ugly
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

} // end Lambert::battin

double
Lambert::bat_SEE(double v)
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
} // end Lambert::bat_SEE

double
Lambert::bat_K(double v)
{
    // d: array (0..20) of Real; -- hardcoded, see astiod.adb
    // XXX ugly
    /*
    double d[20];
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
    d[13]   = 1480.0    / 6075.0;
    d[14]   = 1804.0    / 7047.0;
    d[15]   = 1978.0    / 8091.0;
    d[16]   = 2350.0    / 9207.0;
    d[17]   = 2548.0    / 10395.0;
    d[18]   = 2968.0    / 11655.0;
    d[19]   = 3190.0    / 12987.0;
    d[20]   = 3658.0    / 14391.0;
    */

    const double d[] = {
        0.3333333333333333,
        0.1481481481481481,
        0.2962962962962963,
        0.2222222222222222,
        0.2716049382716049,
        0.2334455667789001,
        0.2641802641802642,
        0.2381766381766382,
        0.2605664488017429,
        0.2407980736154111,
        0.2584238373712058,
        0.2424660685530251,
        0.2436213991769547,
        0.2559954590605932,
        0.2444691632678284,
        0.2552405778212230,
        0.2451178451178451,
        0.2546546546546546,
        0.2456302456302456,
        0.2541866444305469 };

    int     i;
    double  del, delold, term, termold, sum1;

    // process fowards
    sum1 = d[0];
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
} // end Lambert::bat_K
