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
 * $Id: Lambert.cpp,v 1.29 2006/04/06 20:38:08 trs137 Exp $
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

/*
 * UNIVERSAL VARIABLES METHOD
 *
 * Adapted from David Vallado's implementations
 * "Fundamentals of Astrodynamics and Applications"
 */
void
Lambert::universal (const bool L, const int multirev)
{ 
    // Local variables
    const bool longway = L;
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

    /*

    if (true == multirev)
    {
        // For multiple revs only...
        Upper = -0.001 + 16.0 * PI * PI;
        Lower = 0.001 + 4.0 * PI * PI;
    }else{
        // common single revolution case
        Upper = 4.0 * PI * PI;
        Lower = -8.0 * PI;
    }
    */


    if (0 == revs)
    {
        Upper = 4.0 * PI * PI;
        Lower = -8.0 * PI;
    }else{
        Upper = -SMALL + 4.0 * (revs+1)*(revs+1)*PI*PI;
        Lower = SMALL + 4.0 * (revs)*(revs) * PI*PI;
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
            Vo.set3(INF, INF, INF);
            V.set3(INF, INF, INF);
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
        Vo.set3(INF, INF, INF);
        V.set3(INF, INF, INF);
        failure = true;
    } // end if VarA > SMALL

/* Convert from canonical units back to S.I.
    Vo = Vo * ER / TU_SEC;
    V = V * ER / TU_SEC;
*/
}

/*
 * BATTIN'S METHOD
 *
 * Adapted from David Vallado's implementations
 * "Fundamentals of Astrodynamics and Applications"
 */
void
Lambert::battin (void)
{
    // TODO
}

/*
 * Calculate some orbital elements for the xfer arc.
 * a = semimajor axis
 * h = specific angular momentum
 * e = eccentricity
 */
void
Lambert::elements (void)
{
    double VVo;
    VVo = norm(Vo);

    energy = VVo * VVo / 2.0 - 1.0 / norm(Ro);

    // semimajor axis
    a = - 1.0 / (2 * energy);

    // Specific angular momentum
    h = norm(cross(Ro, Vo));

    // Eccentricity  (Canonical MU = 1.0)
    e = sqrt( 1 + (2*energy*h*h));

    // Period of xfer orbit.
   /* 
    if (energy < 0)
    {
        period = 2*PI*sqrt(a*a*a);
    }else{
        period = INF;
    }
   */
}
