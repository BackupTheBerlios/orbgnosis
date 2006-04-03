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
 * $Id: Lambert.cpp,v 1.25 2006/04/03 02:03:28 trs137 Exp $
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

int min_iter = 1000;
int max_iter = 0;
int sum_iter = 0;
int limit = 0;

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
Lambert::universal (void)
{ 
    // Local variables
    const int NumIter = 40;
    int Loops, YNegKtr;
    double VarA, Y, Upper, Lower, CosDeltaNu, F, G, GDot, XOld,
           XOldCubed, PsiOld, PsiNew, C2New, C3New, dtNew;
    double Ro4, R4;
 
    PsiNew = 0.0;
    Vo.toZero();
    V.toZero();

    // Magnitudes of Ro and R
    Ro4 = norm(Ro);
    R4  = norm(R);

    // "Nu" is true anomaly.
    CosDeltaNu = dot(Ro, R) / (Ro4 * R4);

    // We don't care about long way xfers, so...
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

/* Convert from canonical units back to S.I.

    Vo = Vo * ER / TU_SEC;
    V = V * ER / TU_SEC;
*/

    // Count # iterations used for statistics.
    sum_iter = sum_iter + Loops;
    if (max_iter < Loops) max_iter = Loops;
    if (min_iter > Loops) min_iter = Loops;
    if (Loops == NumIter) limit = limit +1;
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
    double energy, VVo;
    VVo = norm(Vo);

    energy = VVo * VVo / 2.0 - 1.0 / norm(Ro);

    // semimajor axis
    a = - 1.0 / (2 * energy);

    // Specific angular momentum
    h = norm(cross(Ro, Vo));

    // Eccentricity  (Canonical MU = 1.0)
    e = sqrt( 1 + (2*energy*h*h));

//    cout << "tof(s)        = " << t*TU_SEC << "\n";
//    cout << "Energy        = " << energy << "\n\n";
//    cout << "Semimajor axis (ER) = " << a << "\n";
//    cout << "h                   = " << h << "\n";
//    cout << "e                   = " << e << "\n";
//    cout << "Perigee (ER)        = " << fabs(a*(1.0+e)) << "\n\n";

cout << t*TU_SEC << ", ";

    
}



int
main(void) {
    cout << "Testing the Lambert solver.\n";
    
    const int problems = 100;
    double t;
    Vector q1, q2;

    Lambert* testcase = new Lambert[problems];

    srand(time(NULL));

    cout << "Generating " << problems << " Lambert's problems\n";
    
/*  TEST 1: Make random problems

    const int prange = 12000;  // -prange to +prange (km)
    const int trange = 5000;   // 0 to +trange (s)
    double a, b, c;
    double    x     = 0.0;     // random double between 0 and 1

   
    for (int i = 0; i < problems; i++)
    {
        cout << ".";
        x = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
        a = (2*x*prange - prange) / ER;
        x = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
        b = (2*x*prange - prange) / ER;
        x = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
        c = (2*x*prange - prange) / ER;

        q.set3(a, b, c);
        testcase[i].setRo(q);

        x = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
        a = (2*x*prange - prange) / ER;
        x = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
        b = (2*x*prange - prange) / ER;
        x = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
        c = (2*x*prange - prange) / ER;

        q.set3(a, b, c);
        testcase[i].setR(q);

        x = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
        t = x*trange / TU_SEC;

        // We're filling thse with CANONICAL UNITS, not S.I.

        testcase[i].sett(t);
    }
*/

/* TEST 2: Solve Curtiss 5.02 for a range of TOFs

    q1.set3(5000.0, 10000.0, 2100.0);
    q1 = q1 / ER;

    q2.set3(-14600.0, 2500.0, 7000.0);
    q2 = q2 / ER;

    double t_max, t_min, t_inc;
    t_min = 60.0 / TU_SEC;  // 1 minute in canonical TU
    // t_max = 1 orbital period for a cirular orbit of radius q1.
    double radius = norm(q1);
    t_max = 2.0 * PI * sqrt(radius*radius*radius);

    // we want (problems) incrments from t_min to t_max.
    t_inc = (t_max-t_min) / (problems-1);

    for (int i = 0; i < problems; i++)
    {
        cout << ".";
        testcase[i].setRo(q1);
        testcase[i].setR(q2);
        testcase[i].sett(t_min+(i*t_inc));
    }

    cout << "\n\nThe problems are ready. Here we go!\n\n";

    for (int i = 0; i < problems; i++)
    {
        // cout << "Problem " << i << ":\n";
        testcase[i].universal();
        testcase[i].elements();
    }

*/

//  TEST 3:
//  r1 = 1.1 (ER);
//  r2 = 1.6 (ER);
//  sweep delta-f (delta-nu) from zero to pi,
//  vary TOF from small to large.

//  Plotting a 2D cartesian of delta-f vs. TOF
//  for delta V should be sort of like a porkchop plot.

    q1.set3(1.1, 0.0, 0.0);  // units of ER
    t = 1461.0 / TU_SEC;

    double f, f_max, f_min, f_inc;
    f_min = SMALL;
    f_max = PI - SMALL;
    f_inc = (f_max-f_min) / (problems-1);

    Vector vc1 (0, sqrt (1/1.1), 0);
    Vector vc2;

    double deltav, t_max, t_min, t_inc;

    t_min = 60.0 / TU_SEC;  // 1 minute in canonical
    // t_max = 1/2 orbital period for a circ radius of 1.6 ER
    t_max = PI * sqrt(1.6*1.6*1.6);
    t_inc = (t_max-t_min) / (problems-1);

    // This will run (problems * problems) times!!!!!

    for (int i = 0; i < problems; i++)
    {
        t = t_min + i*t_inc;

        for (int j = 0; j < problems; j++)
        {
            f = f_min + j*f_inc;

            q2.setX( 1.6*cos(f) );
            q2.setY( 1.6*sin(f) );
            q2.setZ(0.0);

            vc2.setX (-sin(f) * sqrt(1/1.6) );
            vc2.setY (cos(f) * sqrt(1/1.6) );
            vc2.setZ (0.0);


            testcase[j].setRo(q1);
            testcase[j].setR(q2);
            testcase[j].sett(t);

            testcase[j].universal();
            cout << f << ", ";
            testcase[j].elements();

            // delta V 1 in (km/s)
            // cout << norm(testcase[i].getVo() - vc1) * ER / TU_SEC << ", "; 

            // delta V 2 in (km/s)
            // cout << norm(testcase[i].getV() - vc2) * ER / TU_SEC << ", ";

            deltav = ( norm(testcase[j].getVo() - vc1)
                    + norm(testcase[j].getV() - vc2) ) * ER / TU_SEC;

            cout << deltav << "\n";
        }
    }

    delete[] testcase;
    testcase = NULL;

    cout << "\n*****************************************************\n";
    cout << "Tolerance:             " << SMALL << "\n";
    cout << "# Problems:            " << problems*problems << "\n";
    cout << "Max # Iterations:      " << max_iter << "\n";
    cout << "Min # Iterations:      " << min_iter << "\n";
    cout << "Average # Iterations   " << sum_iter/problems/problems << "\n";
    cout << "Exceeded limit  count  " << limit << "\n"; 
    cout << "\n*****************************************************\n";
    return 0;
}
