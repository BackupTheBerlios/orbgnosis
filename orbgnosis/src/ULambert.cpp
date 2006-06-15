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
* $Id: ULambert.cpp,v 1.10 2006/06/15 20:50:33 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*                  David Vallado <valladodl@worldnet.att.net>
*
*/

#include <math.h>
#include "Vec3.h"
#include "Orbgnosis.h"
#include "Stumpff.h"
#include "ULambert.h"
#include <iostream>
using namespace std;

/**
 * The universal variable Lambert constructor with no arguments.
 * All member variables are set to zero.
 */
ULambert::ULambert(void) :
        t(0.0),
        Ro(0.0, 0.0, 0.0),
        R(0.0, 0.0, 0.0),
        Vo(0.0, 0.0, 0.0),
        V(0.0, 0.0, 0.0),
        failure(false)
{
    //cout << "ULambert constructor called with no args.\n";
}

/**
 * The universal variable Lambert constructor with three arguments specified.
 * @param r1in is the initial position.
 * @param r2in is the final position.
 * @param tin is the specified time of flight.
 */
ULambert::ULambert(Vec3 r1in, Vec3 r2in, double tin) :
        t(tin),
        Ro(r1in),
        R(r2in),
        Vo(0.0, 0.0, 0.0),
        V(0.0, 0.0, 0.0),
        failure(false)
{
    // cout << "ULambert constructor called with 2 vectors and 1 time.\n";
}

/**
 * The universal variable Lambert destructor.
 */
ULambert::~ULambert (void)
{
    // cout << "ULambert destructor called \n";
}

/**
 * Sets the initial position vector, Ro.
 */
void
ULambert::setRo (Vec3 vin)
{
    Ro = vin;
}

/**
 * Sets the final position vector, R.
 */
void
ULambert::setR (Vec3 vin)
{
    R = vin;
}

/**
 * Sets the time of flight, t.
 */
void
ULambert::sett (double tin)
{
    t = tin;
}

/**
 * Gets the initial velocity vector, Vo.
 * This is the velocity at the point Ro which satisfies the Lamberts problem.
 */
Vec3
ULambert::getVo (void)
{
    return Vo;
}

/**
 * Gets the final velocity vector, V.
 * This is the velocity at the point R which satisfies the Lamberts problem.
 */
Vec3
ULambert::getV (void)
{
    return V;
}

/**
 * Gets the time of flight, t.
 */
double
ULambert::gett (void)
{
    return t;
}

/**
 * Returns true if failure is true, false if faluire is false.
 * This method is necessary because failure is private.
 */
bool
ULambert::isFailure(void)
{
    if (true == failure)
    {
        return true;
    }
    else
    {
        return false;
    }
}

/**
 * Solves Lamberts Problem using universal variables method.
 * Adapted from David Vallado's Ada implementation in "Fundamentals of
 * Astrodynamics and Applications".
 */
void
ULambert::universal (const bool Lin, const int multirev)
{
    /// True if desired solution is to be the "long way", e.g.
    /// the swept angle is greater than pi or 180 degrees.
    const bool longway = Lin;

    /// The maximum number of iterations allowed.
    const int NumIter = 40;

    int Loops;      //<! loop counter.
    int YNegKtr;    //<! counts how many times Y returned negative.
    double VarA;    //<! see Algorithm #55 of "Fundamentals of Astrodynamics and Applications".
    double Y;       //<! see Algorithm #55.
    double Upper;   //<! upper bound for initial bisection.
    double Lower;   //<! lower bound for initial bisection.
    double CosDeltaNu; //<! the cosine of the swept angle in true anomaly.
    double F;       //<! the f function, see http://scienceworld.wolfram.com/physics/F-andG-Functions.html.
    double G;       //<! the g function, see http://scienceworld.wolfram.com/physics/F-andG-Functions.html.
    double GDot;    //<! the g rate of change of G.
    double XOld;    //<! Greek letter Chi from Algorithm #55.
    double XOldCubed; //<! simply XOld cubed.
    double PsiOld;  //<! Greek letter Psi_n from Algorithm #55
    double PsiNew;  //<! Greek letter Psi_n+1 from Algorithm #55
    double C2New;   //<! Locally stored value of the C2 Stumpff function.
    double C3New;   //<! Locally stored value of the C3 Stumpff function.
    double dtNew;   //<! delta-T from Algoritm #55.


    double Ro4; //<! magnitude of the Vec3 Ro.  The naming convention traces back to Vallado's Ada version which used the 4th element of an array to hold the vector norm of the first 3 elements.
    double R4;  //<! magnitude of the Vec3 R.

    const int revs = multirev; //<! is only roughly related to the number of multiple revolutions when attempting to find multi-rev solutions.

    failure = false;

    PsiNew = 0.0;
    Vo.toZero();
    V.toZero();

    // Magnitudes of Ro and R
    Ro4 = norm(Ro);
    R4 = norm(R);

    // "Nu" is true anomaly.
    CosDeltaNu = dot(Ro, R) / (Ro4 * R4);

    if (true == longway)
    {
        VarA = -sqrt( Ro4 * R4 * (1.0 + CosDeltaNu));
    }
    else
    {
        VarA = sqrt( Ro4 * R4 * (1.0 + CosDeltaNu));
    }

    // Form initial guesses.
    PsiOld = 0.0;
    PsiNew = 0.0;
    XOld = 0.0;
    C2New = 0.5;
    C3New = 1.0 / 6.0;

    // Set up initial bounds for the bisection.
    if (0 == revs)
    {
        Upper = 4.0 * M_PI * M_PI;
        Lower = -8.0 * M_PI;
    }
    else
    {
        Upper = -SMALL + 4.0 * ((revs / 2) + 1) * ((revs / 2) + 1) * M_PI * M_PI;
        Lower = SMALL + 4.0 * (revs / 2) * (revs / 2) * M_PI * M_PI;
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
            }
            else
            {
                Y = Ro4 + R4;
            }
            // Check for negative values of Y.
            if ( (0 < VarA) && (0 > Y))
            {
                YNegKtr = 1;
                while ((0 > Y) && (10 > YNegKtr))
                {
                    PsiNew = 0.8 * (1 / C3New) * (1.0
                                                  - (Ro4 + R4) * sqrt(C2New) / VarA);
                    // Find C2 and C3 functions.
                    C2New = stumpff_C2(PsiNew);
                    C3New = stumpff_C3(PsiNew);
                    PsiOld = PsiNew;
                    Lower = PsiOld;
                    if (fabs(C2New) > SMALL)
                    {
                        Y = Ro4 + R4 - (VarA * (1.0 -
                                                PsiOld * C3New) / sqrt(C2New) );
                    }
                    else
                    {
                        Y = Ro4 + R4;
                    }
                    YNegKtr = YNegKtr + 1;
                } // end while loop.
            } // end if Y neg.

            if (10 > YNegKtr)
            {
                if (fabs(C2New) > SMALL)
                {
                    XOld = sqrt(Y / C2New);
                }
                else
                {
                    XOld = 0.0;
                }
                XOldCubed = XOld * XOld * XOld;
                dtNew = XOldCubed * C3New + VarA * sqrt(Y);

                // Readjust upper and lower bounds.
                if (dtNew < t)
                {
                    Lower = PsiOld;
                }
                else
                {
                    Upper = PsiOld;
                }
                PsiNew = (Upper + Lower) * 0.5;

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
            //cout << "NmIter = " << NumIter << "\n";
            Vo.toInf();
            V.toInf();
            failure = true;
        }
        else
        {

            // Use F and G series to find velocity vectors.
            F = 1.0 - Y / Ro4;
            GDot = 1.0 - Y / R4;
            G = 1.0 / (VarA * sqrt(Y)); // 1 over G

            Vo = (R - F * Ro) * G;
            V = (GDot * R - Ro) * G;
        } // end if answer has converged
    }
    else
    {
        // cout << "Vec3s are 180 degrees apart.\n";
        Vo.toInf();
        V.toInf();
        failure = true;
    } // end if VarA > SMALL

    /* Convert from canonical units back to S.I.
        Vo = Vo * ER / TU_SEC;
        V = V * ER / TU_SEC;
    */
} // end ULambert::universal
