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
* $Id: Tour.cpp,v 1.11 2006/10/02 03:52:53 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*                  Frank Ruskey
*/

#include "Tour.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

/**
 * Tour constructor.  Tell it the # of nodes.
 */

Tour::Tour (int numTargets) try
:
    cols( numTargets + 1 ),
    rows( fact(numTargets) ),
    order( rows ),
    rowCtr( 0 )
{
    // Resize order to have (rows) rows.
    order.resize( rows );
    // Resize each row to hold (cols) columns.

    for ( int i = 0; i < rows; i++ )
        order[i].resize(cols);

    // create filename
    stringstream s;

    string snum;

    s << numTargets;

    s >> snum;

    // e.g. ifstream datafile ("data/5.dat");
    string filename = "data/" + snum + ".dat";

    ifstream datafile (filename.c_str());

    string line;

    if (datafile.is_open())
    {
        for (int r = 0; r < rows; r++)
        {
            getline (datafile, line);

            for (int c = 0; c < cols; c++)
            {
                order[r][c] = atoi(line.substr(c, 1).c_str());
            }
        }

        datafile.close();
    }

    else
        cerr << "Could not open file.\n";
}

catch ( ... )
{
    cerr << "Tour constructor failed.\n";
    exit(1);
}



/**
 * The Tour destructor.
 */
Tour::~Tour ( void )
{
    //cout << "Tour destructor called.\n";
}

void
Tour::printOrder ( void )
{
    for ( int r = 0; r < rows; r++ )
    {
        for ( int c = 0; c < cols; c++ )
        {
            cout << order[r][c] << " ";
        }

        cout << "\n";
    }
}

int
Tour::fact (const int k)
{
    int f = 1;

    for ( int i = 1; i <= k; i++ )
        f = f * i;

    //cout << "The factorial of " << k << " is " << f << ".\n";
    return f;
}

int
Tour::get_target (int r, int c)
{
    return order[r][c];
}
