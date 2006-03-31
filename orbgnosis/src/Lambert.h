/*-
 * Copyright 2006 (c) Ted Stodgell. All rights reserved.
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
 * $Id: Lambert.h,v 1.10 2006/03/31 06:38:47 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 *                  David Vallado <valladodl@worldnet.att.net>
 */

#ifndef _LAMBERT_H_
#define _LAMBERT_H_

class Lambert 
{
    public:

                Lambert (void);         // just zeros

                Lambert (Vector r1in,
                         Vector r2in,
                         double tin);

                virtual ~Lambert (void);

                Lambert (const Lambert&); // copy ctor

                /* Universal Variable method */
                void universal(void);
                /*****************************/

                /* Battin's method */
                // TODO
                /*****************************/

                void setRo (Vector);
                void setR  (Vector);
                void sett  (double);

    private:
        double t;
        Vector Ro;
        Vector R;
        
        Vector Vo;             // initial velocity of xfer arc
        Vector V;              // final velocity of xfer arc

        int NumIter, Loops, YNegKtr;
        double VarA, Y, Upper, Lower, CosDeltaNu, F, G, GDot, XOld,
               XOldCubed, FDot, PsiOld, PsiNew, C2New, C3New, dtNew;

        double Ro4, R4;
};

#endif /* _LAMBERT_H_ */
