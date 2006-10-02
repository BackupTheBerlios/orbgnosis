/* $Id: rank.c,v 1.2 2006/10/02 03:52:53 trs137 Exp $ */
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
/* Rank assignment routine */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to assign rank and crowding distance to a population of size pop_size*/
void assign_rank_and_crowding_distance (population *new_pop)
{
    int flag;
    int i;
    int end;
    int front_size;
    int rank = 1;
    list *orig;
    list *cur;
    list *temp1, *temp2;
    orig = (list *)malloc(sizeof(list));
    cur = (list *)malloc(sizeof(list));
    front_size = 0;
    orig->index = -1;
    orig->parent = NULL;
    orig->child = NULL;
    cur->index = -1;
    cur->parent = NULL;
    cur->child = NULL;
    temp1 = orig;

    for (i = 0; i < popsize; i++)
    {
        insert (temp1, i);
        temp1 = temp1->child;
    }

    do
    {
        if (orig->child->child == NULL)
        {
            new_pop->ind[orig->child->index].rank = rank;
            new_pop->ind[orig->child->index].crowd_dist = INF;
            break;
        }

        temp1 = orig->child;
        insert (cur, temp1->index);
        front_size = 1;
        temp2 = cur->child;
        temp1 = del (temp1);
        temp1 = temp1->child;

        do
        {
            temp2 = cur->child;

            do
            {
                end = 0;
                flag = check_dominance (&(new_pop->ind[temp1->index]), &(new_pop->ind[temp2->index]));

                if (flag == 1)
                {
                    insert (orig, temp2->index);
                    temp2 = del (temp2);
                    front_size--;
                    temp2 = temp2->child;
                }

                if (flag == 0)
                {
                    temp2 = temp2->child;
                }

                if (flag == -1)
                {
                    end = 1;
                }
            }

            while (end != 1 && temp2 != NULL);

            if (flag == 0 || flag == 1)
            {
                insert (cur, temp1->index);
                front_size++;
                temp1 = del (temp1);
            }

            temp1 = temp1->child;
        }

        while (temp1 != NULL);

        temp2 = cur->child;

        do
        {
            new_pop->ind[temp2->index].rank = rank;
            temp2 = temp2->child;
        }

        while (temp2 != NULL);

        assign_crowding_distance_list (new_pop, cur->child, front_size);

        temp2 = cur->child;

        do
        {
            temp2 = del (temp2);
            temp2 = temp2->child;
        }

        while (cur->child != NULL);

        rank += 1;
    }

    while (orig->child != NULL);

    free (orig);

    free (cur);

    return ;
}
