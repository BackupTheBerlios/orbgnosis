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
* $Id: Graph.cpp,v 1.2 2006/10/02 03:52:53 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

#include "Graph.h"
#include "Orbgnosis.h"
#include "Vec3.h"
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

/**
 * Graph constructor.
 * @param n number of targets.
 */
Graph::Graph (int n) :
        node(n),
        numTargets(n)
{
    // Set size for container of trajectories.
    node.resize(numTargets);
    //cout << "Graph constructor called. ";
    //cout << "This graph has " << numTargets << " nodes." << endl;
}

/**
 * Graph copy constuctor.
 */
Graph::Graph (const Graph& copy) :
        node(copy.numTargets),
        numTargets(copy.numTargets)
{
    node.resize(copy.numTargets);

    for (int i = 0; i < numTargets; i++)
        node[i] = copy.node[i];
}

/**
 * Graph destructor
 */
Graph::~Graph (void)
{
    //cout << "Graph destructor called." << endl;
}

/**
 * Prints the entire constellation, one target at a time.
 */
void
Graph::print(void)
{
    cout << "GRAPH PRINTOUT FOR " << numTargets;
    cout << " NODES" << endl;

    for (int i = 0; i < numTargets; i++)
    {
        cout << "Node # " << i << ": " << node[i] << endl;
    }
}

/**
 * Sets every trajectory in the constellation equal to one
 * specified trajectory.
 * @param t is the trajectory to which all trajs in the
 *        constellation will be set.
 */
void
Graph::set_all(Vec3 t)
{
    for (int i = 0; i < numTargets; i++)
    {
        node[i] = t;
    }
}

/**
 * Randomly perturb each traj in a constellation by
 * altering the trajectory's classical elements.
 * @param n specifies the scale of the perturbation.  
 */
void
Graph::noise(double n)
{
    double x, y, z;

    for (int j = 0; j < numTargets; j++)
    {
        x = node[j].getX();
        y = node[j].getY();
        z = node[j].getZ();

        x += 2.0 * n * x * ((double)rand() / ((double)(RAND_MAX) + (double)(1))) - (n * x);
        y += 2.0 * n * y * ((double)rand() / ((double)(RAND_MAX) + (double)(1))) - (n * y);
        z += 2.0 * n * z * ((double)rand() / ((double)(RAND_MAX) + (double)(1))) - (n * z);

        node[j].set3(x, y, z);
    }
}
