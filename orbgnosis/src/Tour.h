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
* $Id: Tour.h,v 1.3 2006/06/15 23:48:13 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

#ifndef _TOUR_H_
#define _TOUR_H_
#include <vector>
using namespace std;

/**
 * A sequence of trajectories that visits each target once.
 */
class Tour
{
    public:
        Tour (int, int);          // constructor
        Tour (const Tour&);       // copy constructor
        Tour& operator = (Tour);  // copy assignment operator
        virtual ~Tour (void);     // destructor

        void printOrder(void);

    private:
        const int cols;            // # of targets, or # of columns
        const int rows;            // # of tours, or # of rows
        vector< vector<int> > order;  // all possible tour orders
        vector<int> p;
        vector<int> pi;
        vector<int> dir;
        int rowCtr;                // row counter
        int temp;

        int factorial(const int);  // returns factorial of an int.

        // Steinhaus-Johnson-Trotter permutation
        void permSJT(int);
        void storeRowSJT(void);
        void exchangeSJT(int, int);

        // Lexicologicl permutation order
        void permLex(void);       // Lexicogical permutation
};

#endif /* _TOUR_H_ */
