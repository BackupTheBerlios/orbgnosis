/* $Id: report.c,v 1.2 2006/10/02 03:52:53 trs137 Exp $ */
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
/* Routines for storing population data into files */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to print the information of a population in a file */
void report_pop (population *pop, FILE *fpt)
{
    int i, j, k;

    for (i = 0; i < popsize; i++)
    {
        for (j = 0; j < nobj; j++)
        {
            fprintf(fpt, "%e\t", pop->ind[i].obj[j]);
        }

        if (ncon != 0)
        {
            for (j = 0; j < ncon; j++)
            {
                fprintf(fpt, "%e\t", pop->ind[i].constr[j]);
            }
        }

        if (nreal != 0)
        {
            for (j = 0; j < nreal; j++)
            {
                fprintf(fpt, "%e\t", pop->ind[i].xreal[j]);
            }
        }

        if (nbin != 0)
        {
            for (j = 0; j < nbin; j++)
            {
                for (k = 0; k < nbits[j]; k++)
                {
                    fprintf(fpt, "%d\t", pop->ind[i].gene[j][k]);
                }
            }
        }

        fprintf(fpt, "%e\t", pop->ind[i].constr_violation);
        fprintf(fpt, "%d\t", pop->ind[i].rank);
        fprintf(fpt, "%e\n", pop->ind[i].crowd_dist);
    }

    return ;
}

/* Function to print the information of feasible and non-dominated population in a file */
void report_feasible (population *pop, FILE *fpt)
{
    int i, j, k;

    for (i = 0; i < popsize; i++)
    {
        if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank == 1)
        {
            for (j = 0; j < nobj; j++)
            {
                fprintf(fpt, "%e\t", pop->ind[i].obj[j]);
            }

            if (ncon != 0)
            {
                for (j = 0; j < ncon; j++)
                {
                    fprintf(fpt, "%e\t", pop->ind[i].constr[j]);
                }
            }

            if (nreal != 0)
            {
                for (j = 0; j < nreal; j++)
                {
                    fprintf(fpt, "%e\t", pop->ind[i].xreal[j]);
                }
            }

            if (nbin != 0)
            {
                for (j = 0; j < nbin; j++)
                {
                    for (k = 0; k < nbits[j]; k++)
                    {
                        fprintf(fpt, "%d\t", pop->ind[i].gene[j][k]);
                    }
                }
            }

            fprintf(fpt, "%e\t", pop->ind[i].constr_violation);
            fprintf(fpt, "%d\t", pop->ind[i].rank);
            fprintf(fpt, "%e\n", pop->ind[i].crowd_dist);
        }
    }

    return ;
}
