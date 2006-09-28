/* $Id: mutation.c,v 1.1 2006/09/28 16:06:27 trs137 Exp $ */
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
/* Mutation routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to perform mutation in a population */
void mutation_pop (population *pop)
{
    int i;

    for (i = 0; i < popsize; i++)
    {
        mutation_ind(&(pop->ind[i]));
    }

    return ;
}

/* Function to perform mutation of an individual */
void mutation_ind (individual *ind)
{
    if (nreal != 0)
    {
        real_mutate_ind(ind);
    }

    if (nbin != 0)
    {
        bin_mutate_ind(ind);
    }

    return ;
}

/* Routine for binary mutation of an individual */
void bin_mutate_ind (individual *ind)
{
    int j, k;
    double prob;

    for (j = 0; j < nbin; j++)
    {
        for (k = 0; k < nbits[j]; k++)
        {
            prob = randomperc();

            if (prob <= pmut_bin)
            {
                if (ind->gene[j][k] == 0)
                {
                    ind->gene[j][k] = 1;
                }

                else
                {
                    ind->gene[j][k] = 0;
                }

                nbinmut += 1;
            }
        }
    }

    return ;
}

/* Routine for real polynomial mutation of an individual */
void real_mutate_ind (individual *ind)
{
    int j;
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;

    for (j = 0; j < nreal; j++)
    {
        if (randomperc() <= pmut_real)
        {
            y = ind->xreal[j];
            yl = min_realvar[j];
            yu = max_realvar[j];
            delta1 = (y - yl) / (yu - yl);
            delta2 = (yu - y) / (yu - yl);
            rnd = randomperc();
            mut_pow = 1.0 / (eta_m + 1.0);

            if (rnd <= 0.5)
            {
                xy = 1.0 - delta1;
                val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (pow(xy, (eta_m + 1.0)));
                deltaq = pow(val, mut_pow) - 1.0;
            }

            else
            {
                xy = 1.0 - delta2;
                val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (pow(xy, (eta_m + 1.0)));
                deltaq = 1.0 - (pow(val, mut_pow));
            }

            y = y + deltaq * (yu - yl);

            if (y < yl)
                y = yl;

            if (y > yu)
                y = yu;

            ind->xreal[j] = y;

            nrealmut += 1;
        }
    }

    return ;
}
