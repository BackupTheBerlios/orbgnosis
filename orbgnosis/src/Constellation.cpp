/*-
* Copyright (c) 2005 Ted Stodgell. All rights reserved.
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
* $Id: Constellation.cpp,v 1.5 2006/09/25 18:44:10 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

#include "Constellation.h"
#include "Orbgnosis.h"
#include "Traj.h"
#include <iostream>
#include <math.h>
#include <vector>
using namespace std;

/**
 * Constellation constructor.
 * @param n number of targets.
 */
Constellation::Constellation (int n) :
    numTargets(n)
{
    // Set size for container of trajectories.
    t10s.resize(numTargets);
    cout << "Constellation constructor called. ";
    cout << "This constellation has " << n << " satellites.";
    cout << endl;
}

/**
 * Constellation copy constuctor.
 */
Constellation::Constellation (const Constellation& copy) :
    numTargets(copy.numTargets)
{
    t10s.resize(copy.numTargets);
    for (int i = 0; i < numTargets; i++)
        t10s[i] = copy.t10s[i];
}

/**
 * Constellation destructor
 */
Constellation::~Constellation (void)
{
    cout << "Constellation destructor called." << endl;
}

/**
 * Prints the entire constellation, one target at a time.
 */
void
Constellation::print(void)
{
    cout << "CONSTELLATION PRINTOUT FOR " << numTargets;
    cout << " SATELLITES" << endl;
    for (int i=0; i<numTargets; i++)
    {
        cout << "Satellite " << i << endl;
        t10s[i].print();
        cout << endl;
    }
}

/**
 * Sets every trajectory in the constellation equal to one
 * specified trajectory.
 * @param t is the trajectory to which all trajs in the
 *        constellation will be set.
 */
void
Constellation::set_all(Traj t)
{
    for (int i=0; i<numTargets; i++)
    {
        t10s[i] = t;
    }
}

/**
 * Given one traj object, fill the rest of the constellation
 * with targets in the same orbit, but distributed evenly around
 * the Earth in terms of mean anomaly.
 * This results in a constellation where the targets are separated
 * by a uniform time interval.
 * @param t is the initial trajectory, t10s[0].
 * All of the other trajectories generated by this function are identical
 * to the inital trajectory except for their anomalies.
 */
void
Constellation::distribute(Traj t)
{
    double M = t.get_M();
    for (int i=0; i < numTargets; i++)
    {
        t10s[i] = t;
        t10s[i].set_M( fmod(M + i*2.0*M_PI/numTargets, 2.0*M_PI));
    }
}

/**
 * Randomly perturb each traj in a constellation by
 * altering the trajectory's classical elements.
 * @param n specifies the scale of the perturbation.  
 */
void
Constellation::noise(double n)
{
    double a, e, i, raan, w, f;
    for (int j=0; j < numTargets; j++)
    {
        a = t10s[j].get_a();
        e = t10s[j].get_e();
        i = t10s[j].get_i();
        raan = t10s[j].get_raan();
        w = t10s[j].get_w();
        f = t10s[j].get_f();

        a += 2.0*n*a*((double)rand()/((double)(RAND_MAX)+(double)(1)))-(n*a);
        e += 2.0*n*e*((double)rand()/((double)(RAND_MAX)+(double)(1)))-(n*e);
        i += 4.0*n*M_PI*((double)rand()/((double)(RAND_MAX)+(double)(1)))-(n*2*M_PI);
        raan += 4.0*n*M_PI*((double)rand()/((double)(RAND_MAX)+(double)(1)))-(n*2*M_PI);
        w += 4.0*n*M_PI*((double)rand()/((double)(RAND_MAX)+(double)(1)))-(n*2*M_PI);
        f += 4.0*n*M_PI*((double)rand()/((double)(RAND_MAX)+(double)(1)))-(n*2*M_PI);

        t10s[j].set_elorb(a, e, i, raan, w, f);
    }
}