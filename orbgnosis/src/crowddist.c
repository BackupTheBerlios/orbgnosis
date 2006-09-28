/* $Id: crowddist.c,v 1.1 2006/09/28 16:06:27 trs137 Exp $ */
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
/* Crowding distance computation routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Routine to compute crowding distance based on ojbective function values when the population in in the form of a list */
void assign_crowding_distance_list (population *pop, list *lst, int front_size)
{
    int **obj_array;
    int *dist;
    int i, j;
    list *temp;
    temp = lst;

    if (front_size == 1)
    {
        pop->ind[lst->index].crowd_dist = INF;
        return ;
    }

    if (front_size == 2)
    {
        pop->ind[lst->index].crowd_dist = INF;
        pop->ind[lst->child->index].crowd_dist = INF;
        return ;
    }

    obj_array = (int **)malloc(nobj * sizeof(int));
    dist = (int *)malloc(front_size * sizeof(int));

    for (i = 0; i < nobj; i++)
    {
        obj_array[i] = (int *)malloc(front_size * sizeof(int));
    }

    for (j = 0; j < front_size; j++)
    {
        dist[j] = temp->index;
        temp = temp->child;
    }

    assign_crowding_distance (pop, dist, obj_array, front_size);
    free (dist);

    for (i = 0; i < nobj; i++)
    {
        free (obj_array[i]);
    }

    free (obj_array);
    return ;
}

/* Routine to compute crowding distance based on objective function values when the population in in the form of an array */
void assign_crowding_distance_indices (population *pop, int c1, int c2)
{
    int **obj_array;
    int *dist;
    int i, j;
    int front_size;
    front_size = c2 - c1 + 1;

    if (front_size == 1)
    {
        pop->ind[c1].crowd_dist = INF;
        return ;
    }

    if (front_size == 2)
    {
        pop->ind[c1].crowd_dist = INF;
        pop->ind[c2].crowd_dist = INF;
        return ;
    }

    obj_array = (int **)malloc(nobj * sizeof(int));
    dist = (int *)malloc(front_size * sizeof(int));

    for (i = 0; i < nobj; i++)
    {
        obj_array[i] = (int *)malloc(front_size * sizeof(int));
    }

    for (j = 0; j < front_size; j++)
    {
        dist[j] = c1++;
    }

    assign_crowding_distance (pop, dist, obj_array, front_size);
    free (dist);

    for (i = 0; i < nobj; i++)
    {
        free (obj_array[i]);
    }

    free (obj_array);
    return ;
}

/* Routine to compute crowding distances */
void assign_crowding_distance (population *pop, int *dist, int **obj_array, int front_size)
{
    int i, j;

    for (i = 0; i < nobj; i++)
    {
        for (j = 0; j < front_size; j++)
        {
            obj_array[i][j] = dist[j];
        }

        quicksort_front_obj (pop, i, obj_array[i], front_size);
    }

    for (j = 0; j < front_size; j++)
    {
        pop->ind[dist[j]].crowd_dist = 0.0;
    }

    for (i = 0; i < nobj; i++)
    {
        pop->ind[obj_array[i][0]].crowd_dist = INF;
    }

    for (i = 0; i < nobj; i++)
    {
        for (j = 1; j < front_size - 1; j++)
        {
            if (pop->ind[obj_array[i][j]].crowd_dist != INF)
            {
                if (pop->ind[obj_array[i][front_size - 1]].obj[i] == pop->ind[obj_array[i][0]].obj[i])
                {
                    pop->ind[obj_array[i][j]].crowd_dist += 0.0;
                }

                else
                {
                    pop->ind[obj_array[i][j]].crowd_dist += (pop->ind[obj_array[i][j + 1]].obj[i] - pop->ind[obj_array[i][j - 1]].obj[i]) / (pop->ind[obj_array[i][front_size - 1]].obj[i] - pop->ind[obj_array[i][0]].obj[i]);
                }
            }
        }
    }

    for (j = 0; j < front_size; j++)
    {
        if (pop->ind[dist[j]].crowd_dist != INF)
        {
            pop->ind[dist[j]].crowd_dist = (pop->ind[dist[j]].crowd_dist) / nobj;
        }
    }

    return ;
}
