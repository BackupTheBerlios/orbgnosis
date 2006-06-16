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
* $Id: Tour.cpp,v 1.5 2006/06/16 20:30:07 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*                  Frank Ruskey
*/

//#include "Orbgnosis.h"
#include "Tour.h"
#include <vector>
#include <iostream>
using namespace std;

/**
 * Tour constructor takes two int args.
 * Generates a matrix of integers rows x columns in size.
 * numTargets is the number of targets.
 * permStyle = 0: Steinhaus-Johnson-Trotter order
 * permStyle = 1: Lexicological order
 */
Tour::Tour (int numTargets, int permStyle) try
    : cols(numTargets),
      rows(factorial(cols)),
      order(rows),
      p(cols+1),
      pi(cols+1),
      dir(cols+1),
      rowCtr(0),
      temp(0)
{
    cout << "Tour constructor called.\n";
    // Resize order to have (rows) rows.
    order.resize(rows);
    for (int i = 0; i < rows; i++)
    {
        // Resize each row to hold (cols) columns.
        order[i].resize(cols);
    }

    switch (permStyle)
    {
        case 0: // Steinhaus-Johnson-Trotter permutations
            // p.resize(cols+1);
            // pi.resize(cols+1);
            // dir.resize(cols+1);
            for (int i = 1; i <= cols; i++)
            {
                dir[i] = -1;
                p[i] = i;
                pi[i] = i;
            }
            permSJT(1);
            break;

        case 1: // Lexicogical permutations
            permLex();
            break;

        default:
            cerr << "Error in Tour constructor: permStyle out of range.\n";
            cerr << "bad permStyle = " << permStyle << ".\n";
            exit(1);
        }
}
catch (...) // std::length_error
{
    cout << "Tour constructor failed.\n";
    exit(1);
}


 
/**
 * The Tour destructor.
 */
Tour::~Tour (void)
{
    cout << "Tour destructor called.\n";
}

/**
 * The Tour copy constructor.
 */
Tour::Tour (const Tour& copy)
    : cols(copy.cols),
      rows(copy.rows),
      order(rows),
      p(cols+1),
      pi(cols+1),
      dir(cols+1),
      rowCtr(0),
      temp(0)
{
    cout << "Tour copy constructor called.\n";
}

/**
 * The Tour copy assignment operator.
 */
Tour&
Tour::operator = (Tour b)
{
    cout << "Tour copy assignment operator was used.\n";
    return *this;
}

/**
 * A private function that returns the factorial of an integer.
 * Only used in Tour member initialization.
 */
int
Tour::factorial (const int k)
{
    int f = 1;
    for (int i = 1; i <= k; i++) f = f * i;
    cout << "The factorial of " << k << " is " << f << ".\n";
    return f;
}

/**
 * The recursive Steinhaus-Johnson-Trotter permutation algorithm.
 * Original algoritm (c) 1995, Frank Ruskey.
 * "can be modified, translated to other languages, etc.
 *  so long as proper acknowledgement is given (author and source)."
 * http://theory.cs.uvic.ca/inf/perm/PermInfo.html
 */
void
Tour::permSJT (int n)
{
    if (n > cols)
    {
        storeRowSJT();
    }
    else
    {
        permSJT (n + 1);
        for (int i = 1; i <= n - 1 ; i++)
        {
            exchangeSJT(n, dir[n]);
            permSJT (n + 1);
        }
        dir[n] = -dir[n];
    }
}

void
Tour::storeRowSJT (void)
{
    for (int i = 0; i < cols; i++)
    {
        order[rowCtr][i] = p[i + 1] - 1; // XXX fixes off-by-one in SJT algo
        //cout << order[rowCtr][i] << " ";
    }
    //cout << "\n";
    rowCtr = rowCtr + 1;
}

void
Tour::exchangeSJT(int x, int d)
{
    temp = p[pi[x] + d];
    p[pi[x]] = temp;
    p[pi[x] + d] = x;
    pi[temp] = pi[x];
    pi[x] = pi[x] + d;
}

void
Tour::printOrder (void)
{
    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            cout << order[r][c] << " ";
        }
        cout << "\n";
    }
}

void
Tour::permLex(void)
{
    // TODO
}
