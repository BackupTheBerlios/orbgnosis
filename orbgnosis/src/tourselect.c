/* $Id: tourselect.c,v 1.2 2006/10/02 03:52:53 trs137 Exp $ */
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
/* Tournamenet Selections routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Routine for tournament selection, it creates a new_pop from old_pop by performing tournament selection and the crossover */
void selection (population *old_pop, population *new_pop)
{
    int *a1, *a2;
    int temp;
    int i;
    int rand;
    individual *parent1, *parent2;
    a1 = (int *)malloc(popsize * sizeof(int));
    a2 = (int *)malloc(popsize * sizeof(int));

    for (i = 0; i < popsize; i++)
    {
        a1[i] = a2[i] = i;
    }

    for (i = 0; i < popsize; i++)
    {
        rand = rnd (i, popsize - 1);
        temp = a1[rand];
        a1[rand] = a1[i];
        a1[i] = temp;
        rand = rnd (i, popsize - 1);
        temp = a2[rand];
        a2[rand] = a2[i];
        a2[i] = temp;
    }

    for (i = 0; i < popsize; i += 4)
    {
        parent1 = tournament (&old_pop->ind[a1[i]], &old_pop->ind[a1[i + 1]]);
        parent2 = tournament (&old_pop->ind[a1[i + 2]], &old_pop->ind[a1[i + 3]]);
        crossover (parent1, parent2, &new_pop->ind[i], &new_pop->ind[i + 1]);
        parent1 = tournament (&old_pop->ind[a2[i]], &old_pop->ind[a2[i + 1]]);
        parent2 = tournament (&old_pop->ind[a2[i + 2]], &old_pop->ind[a2[i + 3]]);
        crossover (parent1, parent2, &new_pop->ind[i + 2], &new_pop->ind[i + 3]);
    }

    free (a1);
    free (a2);
    return ;
}

/* Routine for binary tournament */
individual* tournament (individual *ind1, individual *ind2)
{
    int flag;
    flag = check_dominance (ind1, ind2);

    if (flag == 1)
    {
        return (ind1);
    }

    if (flag == -1)
    {
        return (ind2);
    }

    if (ind1->crowd_dist > ind2->crowd_dist)
    {
        return (ind1);
    }

    if (ind2->crowd_dist > ind1->crowd_dist)
    {
        return (ind2);
    }

    if ((randomperc()) <= 0.5)
    {
        return (ind1);
    }

    else
    {
        return (ind2);
    }
}
