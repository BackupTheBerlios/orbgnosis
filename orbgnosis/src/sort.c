/* $Id: sort.c,v 1.2 2006/10/02 03:52:53 trs137 Exp $ */
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
/* Routines for randomized recursive quick-sort */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Randomized quick sort routine to sort a population based on a particular objective chosen */
void quicksort_front_obj(population *pop, int objcount, int obj_array[], int obj_array_size)
{
    q_sort_front_obj (pop, objcount, obj_array, 0, obj_array_size - 1);
    return ;
}

/* Actual implementation of the randomized quick sort used to sort a population based on a particular objective chosen */
void q_sort_front_obj(population *pop, int objcount, int obj_array[], int left, int right)
{
    int index;
    int temp;
    int i, j;
    double pivot;

    if (left < right)
    {
        index = rnd (left, right);
        temp = obj_array[right];
        obj_array[right] = obj_array[index];
        obj_array[index] = temp;
        pivot = pop->ind[obj_array[right]].obj[objcount];
        i = left - 1;

        for (j = left; j < right; j++)
        {
            if (pop->ind[obj_array[j]].obj[objcount] <= pivot)
            {
                i += 1;
                temp = obj_array[j];
                obj_array[j] = obj_array[i];
                obj_array[i] = temp;
            }
        }

        index = i + 1;
        temp = obj_array[index];
        obj_array[index] = obj_array[right];
        obj_array[right] = temp;
        q_sort_front_obj (pop, objcount, obj_array, left, index - 1);
        q_sort_front_obj (pop, objcount, obj_array, index + 1, right);
    }

    return ;
}

/* Randomized quick sort routine to sort a population based on crowding distance */
void quicksort_dist(population *pop, int *dist, int front_size)
{
    q_sort_dist (pop, dist, 0, front_size - 1);
    return ;
}

/* Actual implementation of the randomized quick sort used to sort a population based on crowding distance */
void q_sort_dist(population *pop, int *dist, int left, int right)
{
    int index;
    int temp;
    int i, j;
    double pivot;

    if (left < right)
    {
        index = rnd (left, right);
        temp = dist[right];
        dist[right] = dist[index];
        dist[index] = temp;
        pivot = pop->ind[dist[right]].crowd_dist;
        i = left - 1;

        for (j = left; j < right; j++)
        {
            if (pop->ind[dist[j]].crowd_dist <= pivot)
            {
                i += 1;
                temp = dist[j];
                dist[j] = dist[i];
                dist[i] = temp;
            }
        }

        index = i + 1;
        temp = dist[index];
        dist[index] = dist[right];
        dist[right] = temp;
        q_sort_dist (pop, dist, left, index - 1);
        q_sort_dist (pop, dist, index + 1, right);
    }

    return ;
}
