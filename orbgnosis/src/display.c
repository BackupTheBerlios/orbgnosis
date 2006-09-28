/* $Id: display.c,v 1.1 2006/09/28 16:06:27 trs137 Exp $ */
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
/* Routines to display the population information using gnuplot */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <unistd.h>

# include "global.h"
# include "rand.h"

/* Function to display the current population for the subsequent generation */
void onthefly_display (population *pop, FILE *gp, int ii)
{
    int i;
    int flag;
    FILE *fpt;
    fpt = fopen("plot.out", "w");
    flag = 0;

    for (i = 0; i < popsize; i++)
    {
        if (pop->ind[i].constr_violation == 0)
        {
            if (choice != 3)
                fprintf(fpt, "%e\t%e\n", pop->ind[i].obj[obj1 - 1], pop->ind[i].obj[obj2 - 1]);
            else
                fprintf(fpt, "%e\t%e\t%e\n", pop->ind[i].obj[obj1 - 1], pop->ind[i].obj[obj2 - 1], pop->ind[i].obj[obj3 - 1]);

            fflush(fpt);

            flag = 1;
        }
    }

    if (flag == 0)
    {
        printf("\n No feasible soln in this pop, hence no display");
    }

    else
    {
        if (choice != 3)
            fprintf(gp, "set title 'Generation #%d'\n unset key\n plot 'plot.out' w points pointtype 6 pointsize 1\n", ii);
        else
            fprintf(gp, "set title 'Generation #%d'\n set view %d,%d\n unset key\n splot 'plot.out' w points pointtype 6 pointsize 1\n", ii, angle1, angle2);

        fflush(gp);
    }

    fclose(fpt);
    return ;
}
