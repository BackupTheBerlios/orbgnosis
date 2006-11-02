/* $Id: crossover.c,v 1.2 2006/11/02 23:12:21 trs137 Exp $ */
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

/* Crossover routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to cross two individuals */
void crossover (individual *parent1, individual *parent2, individual *child1, individual *child2)
{
    if (nreal != 0) realcross (parent1, parent2, child1, child2);
    if (nbin != 0) bincross (parent1, parent2, child1, child2);
    return ;
}

/* Routine for real variable SBX crossover */
void realcross (individual *parent1, individual *parent2, individual *child1, individual *child2)
{
    int i;
    double rand;
    double y1, y2, yl, yu;
    double c1, c2;
    double alpha, beta, betaq;

    if (randomperc() <= pcross_real)
    {
        nrealcross++;

        for (i = 0; i < nreal; i++)
        {
            if (randomperc() <= 0.5 )
            {
                if (fabs(parent1->xreal[i] - parent2->xreal[i]) > EPS)
                {
                    if (parent1->xreal[i] < parent2->xreal[i])
                    {
                        y1 = parent1->xreal[i];
                        y2 = parent2->xreal[i];
                    }

                    else
                    {
                        y1 = parent2->xreal[i];
                        y2 = parent1->xreal[i];
                    }

                    yl = min_realvar[i];
                    yu = max_realvar[i];
                    rand = randomperc();
                    beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1));
                    alpha = 2.0 - pow(beta, -(eta_c + 1.0));

                    if (rand <= (1.0 / alpha))
                    {
                        betaq = pow ((rand * alpha), (1.0 / (eta_c + 1.0)));
                    }

                    else
                    {
                        betaq = pow ((1.0 / (2.0 - rand * alpha)), (1.0 / (eta_c + 1.0)));
                    }

                    c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
                    beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
                    alpha = 2.0 - pow(beta, -(eta_c + 1.0));

                    if (rand <= (1.0 / alpha))
                    {
                        betaq = pow ((rand * alpha), (1.0 / (eta_c + 1.0)));
                    }

                    else
                    {
                        betaq = pow ((1.0 / (2.0 - rand * alpha)), (1.0 / (eta_c + 1.0)));
                    }

                    c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));

                    if (c1 < yl)
                        c1 = yl;

                    if (c2 < yl)
                        c2 = yl;

                    if (c1 > yu)
                        c1 = yu;

                    if (c2 > yu)
                        c2 = yu;

                    if (randomperc() <= 0.5)
                    {
                        child1->xreal[i] = c2;
                        child2->xreal[i] = c1;
                    }

                    else
                    {
                        child1->xreal[i] = c1;
                        child2->xreal[i] = c2;
                    }
                }

                else
                {
                    child1->xreal[i] = parent1->xreal[i];
                    child2->xreal[i] = parent2->xreal[i];
                }
            }

            else
            {
                child1->xreal[i] = parent1->xreal[i];
                child2->xreal[i] = parent2->xreal[i];
            }
        }
    }

    else
    {
        for (i = 0; i < nreal; i++)
        {
            child1->xreal[i] = parent1->xreal[i];
            child2->xreal[i] = parent2->xreal[i];
        }
    }

    return ;
}

/* Routine for two point binary crossover */
void bincross (individual *parent1, individual *parent2, individual *child1, individual *child2)
{
    int i, j;
    double rand;
    int temp, site1, site2;

    for (i = 0; i < nbin; i++)
    {
        rand = randomperc();

        if (rand <= pcross_bin)
        {
            nbincross++;
            site1 = rnd(0, nbits[i] - 1);
            site2 = rnd(0, nbits[i] - 1);

            if (site1 > site2)
            {
                temp = site1;
                site1 = site2;
                site2 = temp;
            }

            for (j = 0; j < site1; j++)
            {
                child1->gene[i][j] = parent1->gene[i][j];
                child2->gene[i][j] = parent2->gene[i][j];
            }

            for (j = site1; j < site2; j++)
            {
                child1->gene[i][j] = parent2->gene[i][j];
                child2->gene[i][j] = parent1->gene[i][j];
            }

            for (j = site2; j < nbits[i]; j++)
            {
                child1->gene[i][j] = parent1->gene[i][j];
                child2->gene[i][j] = parent2->gene[i][j];
            }
        }

        else
        {
            for (j = 0; j < nbits[i]; j++)
            {
                child1->gene[i][j] = parent1->gene[i][j];
                child2->gene[i][j] = parent2->gene[i][j];
            }
        }
    }

    return ;
}
