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
 * $Id: Orbgnosis.cpp,v 1.5 2006/05/05 21:50:40 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 *
 */

#include "Vector.h"
#include "Lambert.h"
#include "Global.h"
#include <math.h>
#include <iostream>
using namespace std;

int
main(void) {
    const int problems = 500;
    double t;
    Vector q1, q2;
    srand(time(NULL));

    Lambert* testcase = new Lambert[problems];

    // Write file for tecplot.

    cout << "TITLE=\"Xfer from 1.1(ER) to 1.2(ER)\"\n";
    cout << "VARIABLES=\"xfer angle (radians)\",\"tof(s)\", \"delta-V(km/s)\"\n";
    cout << "ZONE T=\"MYZONE,BG\", I=" << problems
         << ", J=" << problems << ", F=POINT\n";


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
//  sweep delta-f (delta-nu) from zero to pi,
//  vary TOF from small to large.

//  Plotting a 2D cartesian of delta-f vs. TOF

    q1.set3(1.1, 0.0, 0.0);  // units of ER

    double f, f_max, f_min, f_inc;
    f_min = 2.0 * SMALL;
    f_max = 2.0 * PI + SMALL;
    f_inc = (f_max-f_min) / (problems-1);

    const double r1 = 1.1;
    const double r2 = 1.2;


    Vector vc1 (0, sqrt (1/r1), 0);
    Vector vc2;

    double deltav, short_deltav, long_deltav, t_max, t_min, t_inc;
    bool L;

    t_min = 60.0 / TU_SEC;  // 1 minute in canonical
    t_max = 6.0 * PI * sqrt(r2*r2*r2);
    t_inc = (t_max-t_min) / (problems-1);

    // This will run (problems * problems) times!!!!!

    for (int i = 0; i < problems; i++)
    {
        t = t_min + i*t_inc;

        for (int j = 0; j < problems; j++)
        {
            f = f_min + j*f_inc;

            q2.setX( r2*cos(f) );
            q2.setY( r2*sin(f) );
            q2.setZ(0.0);

            vc2.setX ((-sin(f)) * (sqrt(1/r2)) );
            vc2.setY ((cos(f)) * (sqrt(1/r2)) );
            vc2.setZ (0.0);


            testcase[j].setRo(q1);
            testcase[j].setR(q2);
            testcase[j].sett(t);

            cout << f << ", ";
            cout << t*TU_SEC << ", ";

            // Find minimum delta-V case among many multirev solutions.

            deltav = INF;
            for (int revs = 0; revs < 5; revs++)
            {
                L = false; // short way
                testcase[j].universal(L,revs);
                short_deltav = INF;
                if (!testcase[j].isFailure())
                {
                    short_deltav = ( (norm(testcase[j].getVo() - vc1))
                                  + (norm(testcase[j].getV() - vc2)) )
                                  * ER / TU_SEC;
                }

                L = true; // long way
                testcase[j].universal(L,revs);
                long_deltav = INF;
                if (!testcase[j].isFailure())
                {
                    long_deltav = ( (norm(testcase[j].getVo() - vc1))
                                  + (norm(testcase[j].getV() - vc2)) )
                                  * ER / TU_SEC;
                }

                if (short_deltav < deltav)
                {
                    deltav = short_deltav;
                }

                if (long_deltav < deltav)
                {
                    deltav = long_deltav;
                }
            }
            cout << deltav << "\n";
        }
    }

    delete[] testcase;
    testcase = NULL;

    return 0;
}
