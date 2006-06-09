/*-
 * Copyright 2005 (c) Ted Stodgell. All rights reserved.
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
 * * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * $Id: SJT.cpp,v 1.1 2006/06/09 20:26:11 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 */

#include "SJT.h"
#include <iostream>
using namespace std;

/**
 * The SJT constructor.
 * @param nin is the number of targets in the tour.
 */
SJT::SJT (int nin) : n(nin), m(factorial(n)), currentRow(0)
{
    if ( n > SJT_MAX)
    {
        cerr << "SJT::SJT number of targets exceeds SJT_MAX.\n";
        exit(1);
    }

    try
    {
        p2d = new int* [m];         // allocate double pointer.
        for (int i = 0; i < m; i++)
            p2d[i] = new int [n];   // allocate every i-th row.a
    } catch(...) {
        cerr << "SJT::SJT could not allocate memory for an ";
        cerr << m << " x " << n << " arrray of ints.\n";
        exit(1);
    }

    /* This is just temporary.  Fill array with -999's for no good reason */
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            p2d[i][j] = -999;

    for (int i=1; i<=n; ++i)
    {
        dir[i] = -1; p[i] = i;
        pi[i] = i;
    }
    permutate(1);
}

/**
 * The SJT destructor.
 */
SJT::~SJT (void)
{
    for (int i = 0; i < m; i++)
    {
        delete[] p2d[i];
        p2d[i] = NULL;
    }
    delete[] p2d;
    p2d = NULL;
}

void
SJT::permutate (int a)
{
    int i;
    if (a > n)
    {
        storeRow();
    }else{
        permutate(a+1);
        for (i=1; i <= a-1; ++i)
        {
            exchange(a, dir[a]);
            permutate(a+1);
        }
        dir[a] = -dir[a];
    }
}

void
SJT::storeRow (void)
{
    for (int i = 1; i <=  n; ++i)
    {
        p2d[currentRow][i-1] = p[i]; // SJT algorithm is off by one, this fixes.
    }
    currentRow += 1;
}

void
SJT::exchange (int x, int d)
{
    int z;
    z = p[pi[x]+d];
    p[pi[x]] = z;
    p[pi[x]+d] = x;
    pi[z] = pi[x];
    pi[x] = pi[x]+d;
}

int
SJT::factorial (int n)
{
    int f = 1;
    for (int i=1; i<=n; i++)
        f = f * i;
    return f;
}

void
SJT::print (void)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << p2d[i][j] << " ";
        }
        cout << "\n";
    }
}
