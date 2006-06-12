/*-
* Copyright (c) 2006 Ted Stodgell. All rights reserved.
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
* $Id: SJT.h,v 1.6 2006/06/12 21:22:17 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

#ifndef _TOUR_H_
#define _TOUR_H_

#define SJT_MAX 20

/**
 * All permutations of tour-orders for n targets.
 *
 * The number of permutations is n! (n-factorial).
 * Each tour-order is represented by a key.  For example,
 * the key { 1, 2, 3, 4 } represents one possible tour-order
 * for a system of 4 targets.
 *
 * To ensure efficient operation of the genetic algorithm,
 * successive keys differ only by the exchange of two elements in
 * adjacent positions.  This sequence of permutations is generated
 * by an implementation of the Steinhaus-Johnson-Trotter algorithm.
 * 
 */
class SJT
{
    public:
        SJT (int);              // constructor
        virtual ~SJT (void);            // destructor

        SJT (const SJT&);       // copy constructor
        SJT& operator = (const SJT); // copy assignment operator

        void print (void);         // print the matrix to stdout.

        int getElement(int, int);   // returns p2d[int][int]
        int getRows();              // returns m;
        int getCols();              // returns n;

    private:
        int factorial (int);       // returns factorial
        void permutate (int);       // this is called recursively
        void storeRow (void);      // store 1 tour-order in matrix
        void exchange (int, int);  // exchange 2 elements

        const int n;          //!< The number of targets.
        const int m;          //!< Holds value of n-factorial.a
        int currentRow; //!< counter
        int** p2d;        //!< matrix of ints, m rows by n cols
        int p[SJT_MAX + 1];      //!< a permutation.
        int pi[SJT_MAX + 1];     //!< the permutation's inverse.
        int dir[SJT_MAX + 1];    //!< direction for each element.
        int recursionDepth;

};

#endif /* _SJT_H_ */
