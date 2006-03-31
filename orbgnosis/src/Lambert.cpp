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
 * $Id: Lambert.cpp,v 1.20 2006/03/31 06:38:47 trs137 Exp $
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
}

Lambert::Lambert(Vector r1in, Vector r2in, double tin) :
    t(tin), Ro(r1in), R(r2in) 
{
    // cout << "Lambert constructor called \n";
}

Lambert::~Lambert (void)
{
    // BE SURE TO FREE DYNAMIC STUFF
    // cout << "Lambert destructor called \n";
}

void
Lambert::setRo (Vector vin)
{
    Ro = vin;
    // cout << Ro << "\n";
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


/*
 * UNIVERSAL VARIABLES METHOD
 *
 * Adapted from David Vallado's implementations
 * "Fundamentals of Astrodynamics and Applications"
 */
void
Lambert::universal (void)
{ 
    // Initialize Values
    NumIter = 40;
    PsiNew = 0.0;
    Vo.toZero();
    V.toZero();

    // Magnitudes of Ro and R
    Ro4 = norm(Ro);
    R4  = norm(R);

    // "Nu" is true anomaly in Vallado's notation.
    CosDeltaNu = dot(Ro, R) / (Ro4 * R4);

    // We don't care about long way xfers, so...
    // be careful about t for canonical units.
    // VarA = t * sqrt( Ro4 * R4 * (1.0 + CosDeltaNu));
    VarA = sqrt( Ro4 * R4 * (1.0 + CosDeltaNu));


    // Form initial guesses.
    PsiOld = 0.0;
    PsiNew = 0.0;
    XOld   = 0.0;
    C2New  = 0.5;
    C3New  = 1.0/6.0;

    // Set up initial bounds for the bisection.
    // We don't care about multiple revolutions, so...
    Upper = 4.0 * PI * PI;
    Lower = -8.0 * PI;

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

                //*******************
                // DEBUGGING OUTPUT
/*

                cout << "\nIteration   : " << Loops << "\n";
                cout << "Y(ER)       : " << Y << "\n";
                cout << "Xo(sqrt(ER)): " << XOld << "\n";
                cout << "dtNew(TU)   : " << dtNew << "\n";
                cout << "PsiOld      : " << PsiOld << "\n";
                cout << "PsiNew      : " << PsiNew << "\n";
*/

                //****************************


                // Make sure the first guess isn't too close.
                if (( fabs(dtNew - t) < SMALL) && (1 == Loops))
                {
                    dtNew = t - 1.0;
                }
            } // end if (10 > YNegKtr)
        } // end while loop

        if ( (Loops > NumIter) || (YNegKtr > 10) )
        {
            cout << "Error: Lambert Universal failed to converge. \n";
            if (YNegKtr > 10) cout << "Y is negative\n";
            cout << "NumIter = " << NumIter << "\n";
        }else{

            // Use F and G series to find velocity vectors.
            F    = 1.0 - Y/Ro4;
            GDot = 1.0 - Y/R4;
            G    = 1.0 / (VarA * sqrt(Y)); // 1 over G

            Vo = (R - F*Ro) * G;
            V  = (GDot * R - Ro) * G;
        } // end if answer has converged
    }else{
        cout << "Vectors are 180 degrees apart.\n";
    } // end if VarA > SMALL

//    cout << "\nSolution:\n";
//     cout << "v1(ER/TU) = " << Vo << "\n";
//     cout << "v2(ER/TU) = " << V << "\n";

    Vo = Vo * ER / TU_SEC;
    V = V * ER / TU_SEC;
    cout << "v1(km/s)   = " << Vo << "\n";
    cout << "v2(km/s)   = " << V << "\n";
    cout << "Iterations = " << Loops << "\n\n";

}

    

int
main(void) {
    cout << "Testing the Lambert solver.\n";
    // Vector a(15945.34, 0.0, 0.0);
    // Vector b(12214.83899, 10239.46731, 0.0);
    // double time = 76.0 * 60.0; // minutes

    // Vector a(5000.0, 10000.0, 2100.0);
    // Vector b(-14600.0, 2500.0, 7000.0); 
    // double time = 3600.0; // seconds

    // cout << "double     : " << sizeof(double) << "\n";
    // cout << "long double: " << sizeof(long double) << "\n";

    // cout << "r1(km) = " << a << "\n";
    // cout << "r2(km) = " << b << "\n";
    // cout << "Elapsed time (s) " << time << "\n\n";

    // convert to canonical units
    // a = a / ER;
    // b = b / ER;
    // time = time / TU_SEC;

    // cout << "r1(ER) = " << a << "\n";
    // cout << "r2(ER) = " << b << "\n";
    // cout << "Elapsed time (TU) " << time << "\n\n";

    // Lambert test(a, b, time);

    const int problems = 1000000;

    const int prange = 12000;
    const int trange = 8000;   
    double    x     = 0.0;     // random double between 0 and 1
    double a, b, c, t;
    Vector v;

    Lambert* testcase = new Lambert[problems];

    srand(time(NULL));

    cout << "Generating " << problems << " random Lambert's problems...\n\n";

    for (int i = 0; i < problems; i++)
    {
        x = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
        a = (2*x*prange - prange) / ER;
        x = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
        b = (2*x*prange - prange) / ER;
        x = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
        c = (2*x*prange - prange) / ER;

        v.set3(a, b, c);
        testcase[i].setRo(v);

        x = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
        a = (2*x*prange - prange) / ER;
        x = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
        b = (2*x*prange - prange) / ER;
        x = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
        c = (2*x*prange - prange) / ER;

        v.set3(a, b, c);
        testcase[i].setR(v);

        x = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
        t = (500 + x*trange) / TU_SEC;

        testcase[i].sett(t);
    }

    cout << "The Problems are ready. Here we go....\n";

    for (int i = 0; i < problems; i++)
    {
        cout << "Problem " << i << ":\n";
        testcase[i].universal();
    }

    delete[] testcase;
    testcase = NULL;
        
    return 0;
}
