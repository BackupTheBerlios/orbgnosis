/* $Id: dominance.c,v 1.1 2006/09/28 16:06:27 trs137 Exp $ */
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
/* Domination checking routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Routine for usual non-domination checking
   It will return the following values
   1 if a dominates b
   -1 if b dominates a
   0 if both a and b are non-dominated */

int check_dominance (individual *a, individual *b)
{
    int i;
    int flag1;
    int flag2;
    flag1 = 0;
    flag2 = 0;

    if (a->constr_violation < 0 && b->constr_violation < 0)
    {
        if (a->constr_violation > b->constr_violation)
        {
            return (1);
        }

        else
        {
            if (a->constr_violation < b->constr_violation)
            {
                return ( -1);
            }

            else
            {
                return (0);
            }
        }
    }

    else
    {
        if (a->constr_violation < 0 && b->constr_violation == 0)
        {
            return ( -1);
        }

        else
        {
            if (a->constr_violation == 0 && b->constr_violation < 0)
            {
                return (1);
            }

            else
            {
                for (i = 0; i < nobj; i++)
                {
                    if (a->obj[i] < b->obj[i])
                    {
                        flag1 = 1;

                    }

                    else
                    {
                        if (a->obj[i] > b->obj[i])
                        {
                            flag2 = 1;
                        }
                    }
                }

                if (flag1 == 1 && flag2 == 0)
                {
                    return (1);
                }

                else
                {
                    if (flag1 == 0 && flag2 == 1)
                    {
                        return ( -1);
                    }

                    else
                    {
                        return (0);
                    }
                }
            }
        }
    }
}
