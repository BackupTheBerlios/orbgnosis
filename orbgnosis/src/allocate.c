/* $Id: allocate.c,v 1.1 2006/09/28 16:06:27 trs137 Exp $ */
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
/* Memory allocation and deallocation routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to allocate memory to a population */
void allocate_memory_pop (population *pop, int size)
{
    int i;
    pop->ind = (individual *)malloc(size * sizeof(individual));

    for (i = 0; i < size; i++)
    {
        allocate_memory_ind (&(pop->ind[i]));
    }

    return ;
}

/* Function to allocate memory to an individual */
void allocate_memory_ind (individual *ind)
{
    int j;

    if (nreal != 0)
    {
        ind->xreal = (double *)malloc(nreal * sizeof(double));
    }

    if (nbin != 0)
    {
        ind->xbin = (double *)malloc(nbin * sizeof(double));
        ind->gene = (int **)malloc(nbin * sizeof(int));

        for (j = 0; j < nbin; j++)
        {
            ind->gene[j] = (int *)malloc(nbits[j] * sizeof(int));
        }
    }

    ind->obj = (double *)malloc(nobj * sizeof(double));

    if (ncon != 0)
    {
        ind->constr = (double *)malloc(ncon * sizeof(double));
    }

    return ;
}

/* Function to deallocate memory to a population */
void deallocate_memory_pop (population *pop, int size)
{
    int i;

    for (i = 0; i < size; i++)
    {
        deallocate_memory_ind (&(pop->ind[i]));
    }

    free (pop->ind);
    return ;
}

/* Function to deallocate memory to an individual */
void deallocate_memory_ind (individual *ind)
{
    int j;

    if (nreal != 0)
    {
        free(ind->xreal);
    }

    if (nbin != 0)
    {
        for (j = 0; j < nbin; j++)
        {
            free(ind->gene[j]);
        }

        free(ind->xbin);
        free(ind->gene);
    }

    free(ind->obj);

    if (ncon != 0)
    {
        free(ind->constr);
    }

    return ;
}
