/*- Copyright (c) 2006 Ted Stodgell. All rights reserved.
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
* $Id: Orbgnosis.cpp,v 1.30 2006/10/15 02:31:19 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*
*/

/*
 * Orbgnosis filenames are Capitalized.   (and they're C++)
 * NSGA-II filenames are not capitalized. (and they're ansi C)
 */

#include <iostream>
#include <fstream>
#include <string>
#include "Constellation.h"
#include "Traj.h"
#include "Graph.h"
#include "Kepler.h"
#include "ULambert.h"
#include "Orbgnosis.h"
#include "Tour.h"
#include "Vec3.h"

#include <exception>

//# include <stdio.h>    // from old nsga2r.c
//# include <stdlib.h>   // from old nsga2r.c
# include <math.h>
# include <unistd.h>
# include "global.h"
# include "rand.h"

using namespace std;

// XXX Future work, pick one:
//     1. Wrap all of NSGA-II in a class, or refactor some other way.
//     2. Rewrite NSGA-II from scratch (mandatory for commercial use).
// declare a bunch of global vars for NSGA-II
int nreal;
int nbin;
int nobj;
int ncon;
int popsize;
double pcross_real;
double pcross_bin;
double pmut_real;
double pmut_bin;
double eta_c;
double eta_m;
int ngen;
int nbinmut;
int nrealmut;
int nbincross;
int nrealcross;
int *nbits;
double *min_realvar;
double *max_realvar;
double *min_binvar;
double *max_binvar;
int bitlength;
int choice;
int obj1;
int obj2;
int obj3;
int angle1;
int angle2;


/* OMG THIS ACTUALLY COMPILED */
extern Tour mytour(TARGETS);
extern Graph mygraph(TARGETS + 1);
extern Constellation mycon(TARGETS + 1);  // constellation also has chaser.

/****************************************************************/
/* PROBLEM DEFINITION GOES HERE. DON'T USE problemdef.c         */

// # define wsp1             /* Static wandering salesman problem, 1 objective */
// # define wsp2             /* Static wandering salesman problem, 2 objectives */
# define wsp_astro        /* Dynamic wsp where nodes are satellites */

// Single objective = total length of the tour.
#ifdef wsp1
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    /*
     * Add up the length of the Hamiltonian path, going in the order
     * specified by mytour.
     */
    int start, end; // each edge of the graph has a start node and an end node.
    double d, dtot; // edge weight (distance) and total path distance.
    d = dtot = 0.0;
    /*
     * NOTE! We're using a real coded variable for "key".
     * Real coding the tour permutation key lets us take advantage of the
     * Steinhaus-Johnson-Trotter ordering of the permutation data files.. i.e.
     * adjacent permutation differ by exactly one transposition.  It's a kind
     * of combinatoric gray coding.
     */
    int key;        // the corresponding row number in mytour.
    key = (int)xreal[0];  // convert double to int.

    for (int c = 0; c < mytour.cols - 1; c++) // always start at node #0
    {
        start = mytour.get_target(key, c);    // initially, mytour column 0
        end = mytour.get_target(key, c + 1);  // initially, mytour column 1
        d = norm(mygraph.node[start] - mygraph.node[end]);
        dtot = dtot + d;
    }

    obj[0] = dtot;
    return ;
}

#endif // wsp1

// Two objectives: the horizontal and vertical components of the total tour length.
#ifdef wsp2
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    /*
     * Add up the length of the Hamiltonian path, going in the order
     * specified by mytour.
     */
    int start, end; // each edge of the graph has a start node and an end node.
    double x, xtot; // horizontal component of edge weight (distance) and total path distance.
    double y, ytot; // vertical component edge weight (distance) and total path distance.
    Vec3 edge;
    x = y = xtot = ytot = 0.0;
    /*
     * NOTE! We're using a real coded variable for "key".
     * Real coding the tour permutation key lets us take advantage of the
     * Steinhaus-Johnson-Trotter ordering of the permutation data files.. i.e.
     * adjacent permutation differ by exactly one transposition.  It's a kind
     * of combinatoric gray coding.
     */
    int key;        // the corresponding row number in mytour.
    key = (int)xreal[0];  // convert double to int.

    for (int c = 0; c < mytour.cols - 1; c++) // always start at node #0
    {
        start = mytour.get_target(key, c);    // initially, mytour column 0
        end = mytour.get_target(key, c + 1);  // initially, mytour column 1
        edge = mygraph.node[start] - mygraph.node[end];
        x = fabs(edge.getX());
        y = fabs(edge.getY());
        xtot = xtot + x;
        ytot = ytot + y;
    }

    obj[0] = xtot;
    obj[1] = ytot;
    return ;
}

#endif // wsp2

#ifdef wsp_astro
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    /* CHROMOSOME STRUCTURE:
     * xreal[0] = key  (the tour order)
     * xreal[1] = dwell0
     * xreal[2] = TOF0
     * xreal[3] = dwell1
     * xreal[4] = TOF1
     * xreal[5] = dwell2
     * xreal[6] = TOF2
     * etc...
     *
     * obj[0] = total time of flight
     * obj[1] = total delta-V
     */
    obj[0] = 0.0;
    obj[1] = 0.0;

    int start, end; // each edge of the graph has a start node and an end node.

    /*
     * NOTE! We're using a real coded variable for "key".
     * Real coding the tour permutation key lets us take advantage of the
     * Steinhaus-Johnson-Trotter ordering of the permutation data files.. i.e.
     * adjacent permutation differ by exactly one transposition.  It's a kind
     * of combinatoric gray coding.
     */
    int key;        // the corresponding row number in mytour.
    key = (int)xreal[0];  // convert double to int.
    Vec3 V_start;
    Vec3 V_end;
    Vec3 R_start;
    Vec3 R_end;
    double dv_long, dv_short, dv_best;  // longway and shortway deltaV's
    Traj start_traj;
    Traj end_traj;
    ULambert xfer;

    double* dwell = new double [TARGETS];
    double* TOF = new double [TARGETS];
    double* t_depart = new double [TARGETS];
    double* t_arrive = new double [TARGETS];

    for (int i = 0; i < TARGETS; i++)
    {
        t_depart[i] = 0.0;
        t_arrive[i] = 0.0;
    }

    for (int i = 0; i < TARGETS; i++)
    {
        dwell[i] = xreal[2 * i + 1];
        TOF[i] = xreal[2 * i + 2];
    }

    t_depart[0] = dwell[0];
    t_arrive[0] = dwell[0] + TOF[0];

    for (int i = 1; i < TARGETS; i++)
    {
        t_depart[i] = t_depart[i - 1] + TOF[i - 1] + dwell[i];
        t_arrive[i] = t_arrive[i - 1] + dwell[i] + TOF[i];
    }

    /*
        for (int i = 0; i < TARGETS; i++)
        {
            cout << i;
            cout << ": dwell = " << dwell[i];
            cout << ", TOF = " << TOF[i];
            cout << ", depart = " << t_depart[i];
            cout << ", arrive = " << t_arrive[i];
            cout << endl;
        }
        cout << endl;
    */

    for (int c = 0; c < TARGETS; c++) // always start at node #0
    {
        start = mytour.get_target(key, c);
        end = mytour.get_target(key, c + 1);

        /*
         * Figure out where the chaser is at the time of burn #1 on the current xfer arc.
         */ 
        // mycon.t10s[start] is the trajectory at the beginning of this edge (the c-th edge)
        // t_depart[c] is the time at which we leave upon the c-th transfer arc.
        // start_traj is the location of the start-th node at time t_depart[c].
        // the start-th node is where the chaser is right now.

        try
        {
            start_traj = kepler(mycon.t10s[start], t_depart[c]);
        }

        catch (...)
        {
            cout << "Kepler 1 was bad." << endl;
            exit(1);
        }


        /*
         * Figure out where the current target will be at the time of intercept.. i.e.
         * at the time of burn #2 on the current xfer arc.
         */ 
        // mycon.t10s[end] is the trajectory at the end of this edge.
        // t_arrive[c] is the time at which we will complete the c-th transfer arc.
        // end_traj is the location of the end-th node at time t_arrive[c].
        // the end-th node is where the chaser is headed.

        try
        {
            end_traj = kepler(mycon.t10s[end], t_arrive[c]);
        }

        catch (...)
        {
            cout << "Kepler 2 was bad." << endl;
            exit(1);
        }

        // Extract initial state vector from start_traj.
        R_start = start_traj.get_r();

        V_start = start_traj.get_v();

        // Extract final state vector from end_traj.
        R_end = end_traj.get_r();

        V_end = end_traj.get_v();

        // We already know the TOF for this transfer arc: TOF[c].
        // Get ready for universal Lambert problem.
        xfer.setRo(R_start);

        xfer.setR(R_end);

        xfer.sett(TOF[c]);

        /*
         * Lambert's problem has FOUR solutions:
         * prograde, short-way
         * prograde, long-way
         * retrograde, short-way
         * retrograde, long-way
         *
         * The retrograde solutions are a huge delta-V penalty.
         * Try both long-way and short-way prograde transfers
         * and use whichever is best.
         */ 
        // foo.universal(false, 0 ) means SHORT WAY
        // foo.universal(true, 0 )  means LONG WAY

        dv_best = INF;

        for (int revs = 0; revs < 15; revs++) // multirev kludge
        {
            dv_short = INF;
            dv_long = INF;

            xfer.universal(false, revs);  // short-way

            if (!xfer.isFailure())
            {
                // XXX check hit earth
                dv_short = norm(xfer.getVo() - V_start)
                           + norm(xfer.getV() - V_end);
            }

            xfer.universal(true, revs);  // long-way

            if (!xfer.isFailure())
            {
                // XXX check hit earth
                dv_long = norm(xfer.getVo() - V_start)
                          + norm(xfer.getV() - V_end);
            }

            if (dv_short < dv_best)
                dv_best = dv_short;

            if (dv_long < dv_best)
                dv_best = dv_long;
        }
        // XXX we don't remember which solution (long/short, #revs) was best.
        obj[1] = obj[1] + dv_best;
    }

    // The time-of-flight objective function is quite simple.
    // more general case: obj[0] = t_arrive[TARGETS];
    obj[0] = t_arrive[TARGETS - 1];

    /* CLEAN UP MEMORY */
    delete dwell;

    dwell = NULL;

    delete TOF;

    TOF = NULL;

    delete t_depart;

    t_depart = NULL;

    delete t_arrive;

    t_arrive = NULL;

    // convert units from canonical to km and s.
    obj[0] *= TU_MIN;

    obj[1] *= ERTU;

    // a negative constraint value means a violation.
    // We don't want "Star Trek" style maneuvers, so we will
    // constrain missions that use an obscene amount of delta-V.
    if (obj[1] > 200)  // the cutoff is arbitrary
        constr[0] = -1.0;
    else
        constr[0] = 1.0;

    return ;  // Returning from a void function, just to annoy Brian.
}

#endif // wsp_astro

/****************************************************************/

int main (int argc, char **argv) // arg is a random seed {0...1}
{

    if (argc < 2)
    {
        printf("\n Usage ./orbgnosis random_seed \n");
        exit(1);
    }

    seed = (double)atof(argv[1]);

    if (seed <= 0.0 || seed >= 1.0)
    {
        printf("\n Entered seed value is wrong, seed value must be in (0,1) \n");
        exit(1);
    }

    srand (seed * 2*RAND_MAX);



    Traj mytraj;
    // International Space Station
    //mytraj.set_elorb(1.05354259105, 0.0012287, 0.90124090184, 0.55411411224, 0.46170940032, 1.01);
    //mycon.set_all(mytraj);
    //mycon.distribute(mytraj);
    //for (int i = 0; i <= TARGETS; i++) mycon.t10s[i].set_f(1.0 + i*0.1);
    //mycon.t10s [0].set_f(1.55);
    //mycon.noise(0.05);

    /*
    // 9 sats in 3 planes, spaced widely and evenly
    mycon.t10s[0].set_elorb(1.106, 0.0035, 0.959931, 0.2617993878, 1.5708, 0.1);
    mycon.t10s[1].set_elorb(1.106, 0.0035, 0.959931, 0.2617993878, 1.5708, 0.2617993878);
    mycon.t10s[2].set_elorb(1.106, 0.0035, 0.959931, 0.2617993878, 1.5708, 2.3561944902);
    mycon.t10s[3].set_elorb(1.106, 0.0035, 0.959931, 0.2617993878, 1.5708, 4.4505895926);
    mycon.t10s[4].set_elorb(1.106, 0.0035, 0.959931, 2.3561944902, 1.5708, 0.2617993878);
    mycon.t10s[5].set_elorb(1.106, 0.0035, 0.959931, 2.3561944902, 1.5708, 2.3561944902);
    mycon.t10s[6].set_elorb(1.106, 0.0035, 0.959931, 2.3561944902, 1.5708, 4.4505895926);
    mycon.t10s[7].set_elorb(1.106, 0.0035, 0.959931, 4.4505895926, 1.5708, 0.2617993878);
    mycon.t10s[8].set_elorb(1.106, 0.0035, 0.959931, 4.4505895926, 1.5708, 2.3561944902);
    mycon.t10s[9].set_elorb(1.106, 0.0035, 0.959931, 4.4505895926, 1.5708, 4.4505895926);
    */

    // 9 sats in 3 planes, leader-follower spaced 100km, planes 0.5 deg apart
    mycon.t10s[0].set_elorb(1.106, 0.0035, 1.4835298642, 0.872664626, 1.5708, 0.98);
    mycon.t10s[1].set_elorb(1.106, 0.0035, 1.4835298642, 0.872664626, 1.5708, 1.0);
    mycon.t10s[2].set_elorb(1.106, 0.0035, 1.4835298642, 0.872664626, 1.5708, 1.0156788020);
    mycon.t10s[3].set_elorb(1.106, 0.0035, 1.4835298642, 0.872664626, 1.5708, 1.0313576039);
    //    mycon.t10s[4].set_elorb(1.106, 0.0035, 1.4835298642, 0.0959931089, 1.5708, 1.0);
    //    mycon.t10s[5].set_elorb(1.106, 0.0035, 1.4835298642, 0.0959931089, 1.5708, 1.0156788020);
    //    mycon.t10s[6].set_elorb(1.106, 0.0035, 1.4835298642, 0.0959931089, 1.5708, 1.0313576039);

    //    mycon.t10s[7].set_elorb(1.106, 0.0035, 1.4835298642, 0.1047197551, 1.5708, 1.0);
    //    mycon.t10s[8].set_elorb(1.106, 0.0035, 1.4835298642, 0.1047197551, 1.5708, 1.0156788020);
    //    mycon.t10s[9].set_elorb(1.106, 0.0035, 1.4835298642, 0.1047197551, 1.5708, 1.0313576039);

    mycon.print();
    cout << endl;


    /* WSP2 stuff

    mygraph.set_all(Vec3(100.0, 100.0, 0.0));
    mygraph.noise(1.0);

    //mygraph.print();
    //cout << endl;
    //mytour.printOrder();
    //cout << endl;
    // Do an exhaustive search and find the best wsp2, just to check.
    int start, end;
    Vec3 edge;
    double x, xtot, y, ytot;
    double xbest, ybest;
    int key_xbest, key_ybest;
    x = y = xtot = ytot = 0.0;
    xbest = INF;
    ybest = INF;
    key_xbest = -999;
    key_ybest = -999;
    for (int r = 0; r < mytour.rows; r++)
    {
        xtot = 0.0;
        ytot = 0.0;
        for (int c = 0; c < mytour.cols - 1; c++) // always start at node #0
        {
            start = mytour.get_target(r, c);    // initially, mytour column 0
            end   = mytour.get_target(r, c+1);  // initially, mytour column 1
            edge = mygraph.node[start] - mygraph.node[end];
            x = fabs(edge.getX());
            y = fabs(edge.getY());
            xtot = xtot + x;
            ytot = ytot + y;
        }
        if (xtot < xbest)
        {
            xbest = xtot;
            key_xbest = r;
        }
        if (ytot < ybest)
        {
            ybest = ytot;
            key_ybest = r;
        }
    }
    */

    /*******************************************************/
    int i;
    FILE *fpt1;
    FILE *fpt2;
    FILE *fpt3;
    FILE *fpt4;
    FILE *fpt5;
    FILE *gp;
    population *parent_pop;
    population *child_pop;
    population *mixed_pop;

    fpt1 = fopen("initial_pop.out", "w");
    fpt2 = fopen("final_pop.out", "w");
    fpt3 = fopen("best_pop.out", "w");
    fpt4 = fopen("all_pop.out", "w");
    fpt5 = fopen("params.out", "w");
    fprintf(fpt1, "# This file contains the data of initial population\n");
    fprintf(fpt2, "# This file contains the data of final population\n");
    fprintf(fpt3, "# This file contains the data of final feasible population (if found)\n");
    fprintf(fpt4, "# This file contains the data of all generations\n");
    fprintf(fpt5, "# This file contains information about inputs as read by the program\n");
    printf("\n Enter the problem relevant and algorithm relevant parameters ... ");
    printf("\n Enter the population size (a multiple of 4) : ");
    scanf("%d", &popsize);

    if (popsize < 4 || (popsize % 4) != 0)
    {
        printf("\n population size read is : %d", popsize);
        printf("\n Wrong population size entered, hence exiting \n");
        exit (1);
    }

    printf("\n Enter the number of generations : ");
    scanf("%d", &ngen);

    if (ngen < 1)
    {
        printf("\n number of generations read is : %d", ngen);
        printf("\n Wrong nuber of generations entered, hence exiting \n");
        exit (1);
    }

    printf("\n Enter the number of objectives : ");
    scanf("%d", &nobj);

    if (nobj < 1)
    {
        printf("\n number of objectives entered is : %d", nobj);
        printf("\n Wrong number of objectives entered, hence exiting \n");
        exit (1);
    }

    printf("\n Enter the number of constraints : ");
    scanf("%d", &ncon);

    if (ncon < 0)
    {
        printf("\n number of constraints entered is : %d", ncon);
        printf("\n Wrong number of constraints enetered, hence exiting \n");
        exit (1);
    }

    printf("\n Enter the number of real variables : ");
    scanf("%d", &nreal);

    if (nreal < 0)
    {
        printf("\n number of real variables entered is : %d", nreal);
        printf("\n Wrong number of variables entered, hence exiting \n");
        exit (1);
    }

    if (nreal != 0)
    {
        min_realvar = (double *)malloc(nreal * sizeof(double));
        max_realvar = (double *)malloc(nreal * sizeof(double));

        for (i = 0; i < nreal; i++)
        {
            printf ("\n Enter the lower limit of real variable %d : ", i + 1);
            scanf ("%lf", &min_realvar[i]);
            printf ("\n Enter the upper limit of real variable %d : ", i + 1);
            scanf ("%lf", &max_realvar[i]);

            if (max_realvar[i] <= min_realvar[i])
            {
                printf("\n Wrong limits entered for the min and max bounds of real variable, hence exiting \n");
                exit(1);
            }
        }

        printf ("\n Enter the probability of crossover of real variable (0.6-1.0) : ");
        scanf ("%lf", &pcross_real);

        if (pcross_real < 0.0 || pcross_real > 1.0)
        {
            printf("\n Probability of crossover entered is : %e", pcross_real);
            printf("\n Entered value of probability of crossover of real variables is out of bounds, hence exiting \n");
            exit (1);
        }

        printf ("\n Enter the probablity of mutation of real variables (1/nreal) : ");
        scanf ("%lf", &pmut_real);

        if (pmut_real < 0.0 || pmut_real > 1.0)
        {
            printf("\n Probability of mutation entered is : %e", pmut_real);
            printf("\n Entered value of probability of mutation of real variables is out of bounds, hence exiting \n");
            exit (1);
        }

        printf ("\n Enter the value of distribution index for crossover (5-20): ");
        scanf ("%lf", &eta_c);

        if (eta_c <= 0)
        {
            printf("\n The value entered is : %e", eta_c);
            printf("\n Wrong value of distribution index for crossover entered, hence exiting \n");
            exit (1);
        }

        printf ("\n Enter the value of distribution index for mutation (5-50): ");
        scanf ("%lf", &eta_m);

        if (eta_m <= 0)
        {
            printf("\n The value entered is : %e", eta_m);
            printf("\n Wrong value of distribution index for mutation entered, hence exiting \n");
            exit (1);
        }
    }

    printf("\n Enter the number of binary variables : ");
    scanf("%d", &nbin);

    if (nbin < 0)
    {
        printf ("\n number of binary variables entered is : %d", nbin);
        printf ("\n Wrong number of binary variables entered, hence exiting \n");
        exit(1);
    }

    if (nbin != 0)
    {
        nbits = (int *)malloc(nbin * sizeof(int));
        min_binvar = (double *)malloc(nbin * sizeof(double));
        max_binvar = (double *)malloc(nbin * sizeof(double));

        for (i = 0; i < nbin; i++)
        {
            printf ("\n Enter the number of bits for binary variable %d : ", i + 1);
            scanf ("%d", &nbits[i]);

            if (nbits[i] < 1)
            {
                printf("\n Wrong number of bits for binary variable entered, hence exiting");
                exit(1);
            }

            printf ("\n Enter the lower limit of binary variable %d : ", i + 1);
            scanf ("%lf", &min_binvar[i]);
            printf ("\n Enter the upper limit of binary variable %d : ", i + 1);
            scanf ("%lf", &max_binvar[i]);

            if (max_binvar[i] <= min_binvar[i])
            {
                printf("\n Wrong limits entered for the min and max bounds of binary variable entered, hence exiting \n");
                exit(1);
            }
        }

        printf ("\n Enter the probability of crossover of binary variable (0.6-1.0): ");
        scanf ("%lf", &pcross_bin);

        if (pcross_bin < 0.0 || pcross_bin > 1.0)
        {
            printf("\n Probability of crossover entered is : %e", pcross_bin);
            printf("\n Entered value of probability of crossover of binary variables is out of bounds, hence exiting \n");
            exit (1);
        }

        printf ("\n Enter the probability of mutation of binary variables (1/nbits): ");
        scanf ("%lf", &pmut_bin);

        if (pmut_bin < 0.0 || pmut_bin > 1.0)
        {
            printf("\n Probability of mutation entered is : %e", pmut_bin);
            printf("\n Entered value of probability  of mutation of binary variables is out of bounds, hence exiting \n");
            exit (1);
        }
    }

    if (nreal == 0 && nbin == 0)
    {
        printf("\n Number of real as well as binary variables, both are zero, hence exiting \n");
        exit(1);
    }

    choice = 0;
    printf("\n Do you want to use gnuplot to display the results realtime (0 for NO) (1 for yes) : ");
    scanf("%d", &choice);

    if (choice != 0 && choice != 1)
    {
        printf("\n Entered the wrong choice, hence exiting, choice entered was %d\n", choice);
        exit(1);
    }

    if (choice == 1)
    {
        gp = popen(GNUPLOT_COMMAND, "w");

        if (gp == NULL)
        {
            printf("\n Could not open a pipe to gnuplot, check the definition of GNUPLOT_COMMAND in file global.h\n");
            printf("\n Edit the string to suit your system configuration and rerun the program\n");
            exit(1);
        }

        if (nobj == 2)
        {
            printf("\n Enter the objective for X axis display : ");
            scanf("%d", &obj1);

            if (obj1 < 1 || obj1 > nobj)
            {
                printf("\n Wrong value of X objective entered, value entered was %d\n", obj1);
                exit(1);
            }

            printf("\n Enter the objective for Y axis display : ");
            scanf("%d", &obj2);

            if (obj2 < 1 || obj2 > nobj)
            {
                printf("\n Wrong value of Y objective entered, value entered was %d\n", obj2);
                exit(1);
            }

            obj3 = -1;
        }

        else
        {
            printf("\n #obj > 2, 2D display or a 3D display ?, enter 2 for 2D and 3 for 3D :");
            scanf("%d", &choice);

            if (choice != 2 && choice != 3)
            {
                printf("\n Entered the wrong choice, hence exiting, choice entered was %d\n", choice);
                exit(1);
            }

            if (choice == 2)
            {
                printf("\n Enter the objective for X axis display : ");
                scanf("%d", &obj1);

                if (obj1 < 1 || obj1 > nobj)
                {
                    printf("\n Wrong value of X objective entered, value entered was %d\n", obj1);
                    exit(1);
                }

                printf("\n Enter the objective for Y axis display : ");
                scanf("%d", &obj2);

                if (obj2 < 1 || obj2 > nobj)
                {
                    printf("\n Wrong value of Y objective entered, value entered was %d\n", obj2);
                    exit(1);
                }

                obj3 = -1;
            }

            else
            {
                printf("\n Enter the objective for X axis display : ");
                scanf("%d", &obj1);

                if (obj1 < 1 || obj1 > nobj)
                {
                    printf("\n Wrong value of X objective entered, value entered was %d\n", obj1);
                    exit(1);
                }

                printf("\n Enter the objective for Y axis display : ");
                scanf("%d", &obj2);

                if (obj2 < 1 || obj2 > nobj)
                {
                    printf("\n Wrong value of Y objective entered, value entered was %d\n", obj2);
                    exit(1);
                }

                printf("\n Enter the objective for Z axis display : ");
                scanf("%d", &obj3);

                if (obj3 < 1 || obj3 > nobj)
                {
                    printf("\n Wrong value of Z objective entered, value entered was %d\n", obj3);
                    exit(1);
                }

                printf("\n You have chosen 3D display, hence location of eye required \n");
                printf("\n Enter the first angle (an integer in the range 0-180) (if not known, enter 60) :");
                scanf("%d", &angle1);

                if (angle1 < 0 || angle1 > 180)
                {
                    printf("\n Wrong value for first angle entered, hence exiting \n");
                    exit(1);
                }

                printf("\n Enter the second angle (an integer in the range 0-360) (if not known, enter 30) :");
                scanf("%d", &angle2);

                if (angle2 < 0 || angle2 > 360)
                {
                    printf("\n Wrong value for second angle entered, hence exiting \n");
                    exit(1);
                }
            }
        }
    }

    printf("\n Input data successfully entered, now performing initialization \n");
    fprintf(fpt5, "\n Population size = %d", popsize);
    fprintf(fpt5, "\n Number of generations = %d", ngen);
    fprintf(fpt5, "\n Number of objective functions = %d", nobj);
    fprintf(fpt5, "\n Number of constraints = %d", ncon);
    fprintf(fpt5, "\n Number of real variables = %d", nreal);

    if (nreal != 0)
    {
        for (i = 0; i < nreal; i++)
        {
            fprintf(fpt5, "\n Lower limit of real variable %d = %e", i + 1, min_realvar[i]);
            fprintf(fpt5, "\n Upper limit of real variable %d = %e", i + 1, max_realvar[i]);
        }

        fprintf(fpt5, "\n Probability of crossover of real variable = %e", pcross_real);
        fprintf(fpt5, "\n Probability of mutation of real variable = %e", pmut_real);
        fprintf(fpt5, "\n Distribution index for crossover = %e", eta_c);
        fprintf(fpt5, "\n Distribution index for mutation = %e", eta_m);
    }

    fprintf(fpt5, "\n Number of binary variables = %d", nbin);

    if (nbin != 0)
    {
        for (i = 0; i < nbin; i++)
        {
            fprintf(fpt5, "\n Number of bits for binary variable %d = %d", i + 1, nbits[i]);
            fprintf(fpt5, "\n Lower limit of binary variable %d = %e", i + 1, min_binvar[i]);
            fprintf(fpt5, "\n Upper limit of binary variable %d = %e", i + 1, max_binvar[i]);
        }

        fprintf(fpt5, "\n Probability of crossover of binary variable = %e", pcross_bin);
        fprintf(fpt5, "\n Probability of mutation of binary variable = %e", pmut_bin);
    }

    fprintf(fpt5, "\n Seed for random number generator = %e", seed);
    bitlength = 0;

    if (nbin != 0)
    {
        for (i = 0; i < nbin; i++)
        {
            bitlength += nbits[i];
        }
    }

    fprintf(fpt1, "# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n", nobj, ncon, nreal, bitlength);
    fprintf(fpt2, "# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n", nobj, ncon, nreal, bitlength);
    fprintf(fpt3, "# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n", nobj, ncon, nreal, bitlength);
    fprintf(fpt4, "# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n", nobj, ncon, nreal, bitlength);
    nbinmut = 0;
    nrealmut = 0;
    nbincross = 0;
    nrealcross = 0;
    parent_pop = (population *)malloc(sizeof(population));
    child_pop = (population *)malloc(sizeof(population));
    mixed_pop = (population *)malloc(sizeof(population));
    allocate_memory_pop (parent_pop, popsize);
    allocate_memory_pop (child_pop, popsize);
    allocate_memory_pop (mixed_pop, 2*popsize);
    randomize();
    initialize_pop (parent_pop);
    printf("\n Initialization done, now performing first generation");
    decode_pop(parent_pop);
    evaluate_pop (parent_pop);
    assign_rank_and_crowding_distance (parent_pop);
    report_pop (parent_pop, fpt1);
    fprintf(fpt4, "# gen = 1\n");
    report_pop(parent_pop, fpt4);
    printf("\n gen = 1");
    fflush(stdout);

    if (choice != 0)
        onthefly_display (parent_pop, gp, 1);

    fflush(fpt1);

    fflush(fpt2);

    fflush(fpt3);

    fflush(fpt4);

    fflush(fpt5);

    sleep(1);   /* XXX why? */

    for (i = 2; i <= ngen; i++)
    {
        selection (parent_pop, child_pop);
        mutation_pop (child_pop);
        decode_pop(child_pop);
        evaluate_pop(child_pop);
        merge (parent_pop, child_pop, mixed_pop);
        fill_nondominated_sort (mixed_pop, parent_pop);
        /* Comment following four lines if information for all
        generations is not desired, it will speed up the execution */
        fprintf(fpt4, "# gen = %d\n", i);
        report_pop(parent_pop, fpt4);
        fflush(fpt4);

        if (choice != 0)
            onthefly_display (parent_pop, gp, i);

        /* Comment the four lines above for no display */
        printf("\n gen = %d", i);

        // sleep(1);
    }

    printf("\n Generations finished, now reporting solutions");
    report_pop(parent_pop, fpt2);
    report_feasible(parent_pop, fpt3);

    if (nreal != 0)
    {
        fprintf(fpt5, "\n Number of crossover of real variable = %d", nrealcross);
        fprintf(fpt5, "\n Number of mutation of real variable = %d", nrealmut);
    }

    if (nbin != 0)
    {
        fprintf(fpt5, "\n Number of crossover of binary variable = %d", nbincross);
        fprintf(fpt5, "\n Number of mutation of binary variable = %d", nbinmut);
    }

    fflush(stdout);
    fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    fflush(fpt4);
    fflush(fpt5);
    fclose(fpt1);
    fclose(fpt2);
    fclose(fpt3);
    fclose(fpt4);
    fclose(fpt5);

    if (choice != 0)
    {
        pclose(gp);
    }

    if (nreal != 0)
    {
        free (min_realvar);
        free (max_realvar);
    }

    if (nbin != 0)
    {
        free (min_binvar);
        free (max_binvar);
        free (nbits);
    }

    deallocate_memory_pop (parent_pop, popsize);
    deallocate_memory_pop (child_pop, popsize);
    deallocate_memory_pop (mixed_pop, 2*popsize);
    free (parent_pop);
    free (child_pop);
    free (mixed_pop);
    printf("\n Routine successfully exited \n");

    /*
        cout << "The shortest horizontal tour was # " << key_xbest << ", length = ";
        cout << xbest << endl;
        cout << "The shortest vertical tour was # " << key_ybest << ", length = ";
        cout << ybest << endl;
    */


    return EXIT_SUCCESS;
}
