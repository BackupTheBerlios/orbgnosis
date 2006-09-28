/* $Id: rand.c,v 1.1 2006/09/28 16:06:27 trs137 Exp $ */
/*************************************************************************
 * Copyright Notice:                                                     *
 * Source code for random number generator (files rand.h & rand.c) has   *
 * been taken from : sga.c (C) David E. Goldberg 1986.                   *
 * Entire source code (other than mentioned above) present in this       *
 * directory has been developed (from scratch) at Kanpur Genetic         *
 * Algorithms Laboratory (KanGAL, IIT Kanpur) and is the property of its *
 * authors. (C) Dr. Kalyanmoy Deb 2005.                                  *
 *                                                                       *
 * Disclaimer Notice:                                                    *
 * These codes have been developed for research purpose and are in a     *
 * process of continous change, and therefore bugs may exist.            *
 * In any case, developers of the codes do not take any responsibility   *
 * of any malfunction, although they have been tested on many test       *
 * problems. These codes have been tested on Mandrake and Gentoo linux   *
 * (kernel version 2.6.x and gcc version 3.3.x). Any bug or error may    *
 * kindly be communicated to deb@iitk.ac.in. Commercial use of these     *
 * codes is strictly prohibited without the knowledge of developers. For *
 * academic use, these can be used or modified  at will, however an      *
 * acknowledgement of developers at appropriate places would be highly   *
 * appreciated.                                                          *
 *************************************************************************/
/* Definition of random number generation routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

double seed;
double oldrand[55];
int jrand;

/* Get seed number for random and start it up */
void randomize()
{
    int j1;

    for (j1 = 0; j1 <= 54; j1++)
    {
        oldrand[j1] = 0.0;
    }

    jrand = 0;
    warmup_random (seed);
    return ;
}

/* Get randomize off and running */
void warmup_random (double seed)
{
    int j1, ii;
    double new_random, prev_random;
    oldrand[54] = seed;
    new_random = 0.000000001;
    prev_random = seed;

    for (j1 = 1; j1 <= 54; j1++)
    {
        ii = (21 * j1) % 54;
        oldrand[ii] = new_random;
        new_random = prev_random - new_random;

        if (new_random < 0.0)
        {
            new_random += 1.0;
        }

        prev_random = oldrand[ii];
    }

    advance_random ();
    advance_random ();
    advance_random ();
    jrand = 0;
    return ;
}

/* Create next batch of 55 random numbers */
void advance_random ()
{
    int j1;
    double new_random;

    for (j1 = 0; j1 < 24; j1++)
    {
        new_random = oldrand[j1] - oldrand[j1 + 31];

        if (new_random < 0.0)
        {
            new_random = new_random + 1.0;
        }

        oldrand[j1] = new_random;
    }

    for (j1 = 24; j1 < 55; j1++)
    {
        new_random = oldrand[j1] - oldrand[j1 - 24];

        if (new_random < 0.0)
        {
            new_random = new_random + 1.0;
        }

        oldrand[j1] = new_random;
    }
}

/* Fetch a single random number between 0.0 and 1.0 */
double randomperc()
{
    jrand++;

    if (jrand >= 55)
    {
        jrand = 1;
        advance_random();
    }

    return ((double)oldrand[jrand]);
}

/* Fetch a single random integer between low and high including the bounds */
int rnd (int low, int high)
{
    int res;

    if (low >= high)
    {
        res = low;
    }

    else
    {
        res = low + (randomperc() * (high - low + 1));

        if (res > high)
        {
            res = high;
        }
    }

    return (res);
}

/* Fetch a single random real number between low and high including the bounds */
double rndreal (double low, double high)
{
    return (low + (high - low)*randomperc());
}
