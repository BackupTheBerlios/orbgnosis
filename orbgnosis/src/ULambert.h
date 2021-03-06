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
* $Id: ULambert.h,v 1.10 2006/10/02 03:52:53 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*                  David Vallado <valladodl@worldnet.att.net>
*/

#ifndef _ULAMBERT_H_
#define _ULAMBERT_H_

/**
 * The universal variables method of solving the Lambert's problem.
 * Solving Lamberts problem with this technique is difficult because the
 * iteration isn't always well-behaved.  Simple Newtonian iteration may fail
 * to converte on diffucult hyperbolic trajectories.  Therefore, this class
 * uses a simple bisection technique which works well on all trajectory types
 * and is only slightly slower.  The upper and lower bounds are set accordingly
 * so that the solution converges to the desired solution, but only the
 * zero-revolution case(s) are garuanteed.
 */

class ULambert
{

    public:

        ULambert ( void );         // just zeros

        ULambert ( Vec3 r1in,        //!< initial position.
                   Vec3 r2in,        //!< final position
                   double tin ); //!< time of flight.

        virtual ~ULambert ( void );

        // Universal Variable method
        void universal( const bool, const int );

        void setRo ( Vec3 );
        void setR ( Vec3 );
        void sett ( double );

        Vec3 getVo ( void );
        Vec3 getV ( void );
        double gett ( void );

        bool isFailure( void );

    private:

        // Inputs:
        double t;   //!< specified time of flight from Ro to R.
        Vec3 Ro;  //!< initial position vector.
        Vec3 R;   //!< final position vector.

        // Results:
        Vec3 Vo;  //!< initial velocity of at start of the transfer arc.
        Vec3 V;   //!< final velocity at the end of transfer arc.

        bool failure;   //!< is true if the solution fails to converge.
};

#endif /* _ULAMBERT_H_ */
