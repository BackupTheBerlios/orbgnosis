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
* * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
* SUCH DAMAGE.
*
* $Id: SJT.cpp,v 1.7 2006/06/13 23:19:18 trs137 Exp $
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
SJT::SJT (int nin) try
:
    n(nin),
    m(factorial(n)),
    currentRow(0),
    p2d(new int*[m])
{
    for (int i = 0; i < m; i++)
        p2d[i] = new int [n];   // allocate every i-th row.a

    for (int i = 1; i <= n; ++i)
    {
        dir[i] = -1;
        p[i] = i;
        pi[i] = i;
    }
    permutate(1);
}
catch (...)
{
    std::cerr << "Couldn't create SJT.\n";
    exit(1);
}

/**
 * The SJT destructor.
 */
SJT::~SJT (void)
{
    try
    {
        for (int i = 0; i < m; i++)
        {
            delete[] p2d[i];
            p2d[i] = NULL;
        }
        delete[] p2d;
        p2d = NULL;
    }
    catch (...)
    {
        std::cerr << "SJT::~SJT could not free and delete memory.\n";
        exit(1);
    }
}

void
SJT::permutate (const int &a)
{
    if (a > n)
    {
        storeRow();
    }
    else
    {
        permutate(a + 1);
        for (int i = 1; i <= a - 1; ++i)
        {
            exchange(a, dir[a]);
            permutate(a + 1);
        }
        dir[a] = -dir[a];
    }
}

void
SJT::storeRow (void)
{
    for (int i = 0; i < n; ++i)
    {
        p2d[currentRow][i] = p[i + 1] - 1; // XXX SJT algorithm is off by one, this fixes.
    }
    currentRow += 1;
}

void
SJT::exchange (const int &x, int &d)
{
    int z;
    z = p[pi[x] + d];
    p[pi[x]] = z;
    p[pi[x] + d] = x;
    pi[z] = pi[x];
    pi[x] = pi[x] + d;
}

int
SJT::factorial (const int &x)
{
    if ((SJT_MAX < x) || ( 1 > x ))
    {
        std::cerr << "error, SJT::factorial tried to find factorial of " << x;
        std::cerr << "\n but is restricted to positive integers of ";
        std::cerr << SJT_MAX << " or less.\n";
        exit(1);
    }
    int f = 1;
    for (int i = 1; i <= x; i++)
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
            std::cout << p2d[i][j] << " ";
        }
        std::cout << "\n";
    }
}

int
SJT::getElement (int row, int col)
{
    return p2d[row][col];
}

int
SJT::getRows (void)
{
    return m;
}

int
SJT::getCols (void)
{
    return n;
}
