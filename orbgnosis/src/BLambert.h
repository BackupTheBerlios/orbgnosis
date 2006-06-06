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
 * $Id: BLambert.h,v 1.4 2006/06/06 01:57:40 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 *                  David Vallado <valladodl@worldnet.att.net>
 */

#ifndef _BLAMBERT_H_
#define _BLAMBERT_H_

/**
 * BLambert is Battin's method of solving the Lambert's problem.
 * It is robust and does not suffer from from convergence failure.
 * This implementation is adapted from Vallado's Ada procedure.
 */
class BLambert 
{
    public:

            BLambert (void);         // just zeros

            BLambert (Vector r1in,
                         Vector r2in,
                         double tin);

            virtual ~BLambert (void);

            BLambert (const BLambert&); // copy ctor

            // Battin's method
            void    battin(void);
            double  bat_SEE(double);
            double  bat_K(double);

            void setRo (Vector);
            void setR  (Vector);
            void sett  (double);

            Vector getVo (void);
            Vector getV  (void);
            double gett  (void);

            bool isFailure(void);

    private:

        // Inputs:
        double t;   //<! specified time of flight from Ro to R.
        Vector Ro;  //<! initial position vector.
        Vector R;   //<! final position vector.

        // Results:        
        Vector Vo;  //<! initial velocity of at start of the transfer arc.
        Vector V;   //<! final velocity at the end of transfer arc.

        bool failure;   //<! is true if the solution fails to converge.
};

#endif /* _BLAMBERT_H_ */
