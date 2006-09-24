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
* $Id: Constellation.h,v 1.2 2006/09/24 23:57:51 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

#ifndef _CONSTELLATION_H_
#define _CONSTELLATION_H_
#include "Traj.h" 
#include <vector>
using namespace std;

/**
 * A group of target satellites to be visited.
 */
class Constellation
{
    public:
        Constellation (int);                    // constructor
        Constellation (const Constellation&);   // copy ctor

        virtual ~Constellation (void);          // destructor

        void print (void);   // Prints the entire constellation..


        // Methods for generating realistic constellations...
        void set_all(Traj);    // Sets every traj exactly the same.
        void distribute(Traj); // Spreads targets evenly around in one orbit.
        void noise(double);    // Perturbs each target slightly

        // Recommended procedure to make a constellation:
        // 1. Pick a characteristic element set with any arbitrary true anomaly.
        // 2. Use distribute() to populate the plane, using that elset.
        // 3. If desired, use Traj::set_raan() and/or Traj::set_w to crank some
        //    of the orbits in to different plane(s).
        // 4. If desired, use noise() to "jiggle" each element by a small scaling
        //    factor.

        vector<Traj> t10s;          // t10s is short for "trajectories"
        const int numTargets;       // # of satellites in constellation
};

#endif /* _CONSTELLATION_H_ */
