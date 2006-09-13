/*-
* Copyright (c) 2005 Ted Stodgell. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
* OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
* SUCH DAMAGE.
*
* $Id: problemdef.h,v 1.7 2006/09/13 02:01:15 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

/**
 * @file
 * This header file defines the test_problem() function used by NSGA-2.
 * Please see NSGA-2 docs for more examples.
 */

#ifndef _PROBLEMDEF_H_
#define _PROBLEMDEF_H_
#include "Orbgnosis.h"

/**
 * Interface with NSGA-2
 */
void test_problem ( double *xreal, double *xbin, int **gene, double *obj, double *constr )
{
    obj[ 0 ] =      /* TODO */
        obj[ 1 ] =      /* TODO */
            obj[ 2 ] =      /* TODO */
                constr[ 0 ] =      /* TODO */
                    constr[ 1 ] =      /* TODO */
                        constr[ 2 ] =      /* TODO */
                            constr[ 3 ] =      /* TODO */
                                return ;
}

#endif /* _PROBLEMDEF_H_ */
