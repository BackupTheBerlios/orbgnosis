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
 * $Id: Lambert.h,v 1.15 2006/04/06 20:38:08 trs137 Exp $
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

                // Universal Variable method
                void universal(const bool, const int);

                // Battin's method
                void battin(void);

                // Calculate elements from (Vo, V)
                void elements(void);

                void setRo (Vector);
                void setR  (Vector);
                void sett  (double);
                //void setbad(void);

                Vector getVo (void);
                Vector getV  (void);
                double gett  (void);

                bool failure;
                double energy, e;

    private:

        // Inputs:
        double t;
        Vector Ro;
        Vector R;

        // Results:        
        Vector Vo;             // initial velocity of xfer arc
        Vector V;              // final velocity of xfer arc

        // TODO: orbital elements for xfer arc
        double h, a,;
};

#endif /* _LAMBERT_H_ */
