/* $Id: decode.c,v 1.1 2006/09/28 16:06:27 trs137 Exp $ */
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
/* Routines to decode the population */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to decode a population to find out the binary variable values based on its bit pattern */
void decode_pop (population *pop)
{
    int i;

    if (nbin != 0)
    {
        for (i = 0; i < popsize; i++)
        {
            decode_ind (&(pop->ind[i]));
        }
    }

    return ;
}

/* Function to decode an individual to find out the binary variable values based on its bit pattern */
void decode_ind (individual *ind)
{
    int j, k;
    int sum;

    if (nbin != 0)
    {
        for (j = 0; j < nbin; j++)
        {
            sum = 0;

            for (k = 0; k < nbits[j]; k++)
            {
                if (ind->gene[j][k] == 1)
                {
                    sum += pow(2, nbits[j] - 1 - k);
                }
            }

            ind->xbin[j] = min_binvar[j] + (double)sum * (max_binvar[j] - min_binvar[j]) / (double)(pow(2, nbits[j]) - 1);
        }
    }

    return ;
}
