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
* $Id: HitEarth.h,v 1.1 2006/10/15 04:12:13 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*                  David Vallado <valldodl@worldnet.att.net>
*/

#include "Orbgnosis.h"
// maybe... #include "Traj.h"
#include "Vec3.h"
#include <iostream>
#include <math.h>
using namespace std;

/**
 * Returns true if the given trajectory segment passes inside the Earth.
 * @param r1 initial position, type Vec3.
 * @param r2 final position, type Vec3.
 * @param v1 initial velocity, type Vec3.
 * @param v2 final velocity, type Vec3.
 */
inline bool
hit_Earth (const Vec3 r1, const Vec3 r2, const Vec3 v1, const Vec3 v2)
{
    cout << "hit_Earth() is checking...";
    // Are the inital or final points inside the Earth?
    if ((norm(r1) < 1.0) || (norm(r2) < 1.0))
    {
        cout << "HIT!" << endl;
        return true;
    }
    // From Vallado's book: Recall the properties of the flight path angle
    // for an orbit.  The FPA is positive as the satellite travels from perigee
    // to apogee and negative on the return.  The dot product of the position
    // and velocity vectors gives this change of sign without the trigonometric
    // calculations required to calculate the FPA.
    // ... (snip) ... The only case we need to check is when the initial dot
    // product is negative and the final dot product is positive, indicating
    // that perigee occurred during the transfer arc.

    if ((dot(r1, v1) < 0.0) && (dot(r2, v2) > 0.0))
    {
        double h;           // angular momentum of arc
        double ksi;         // specific mechanical energy of arc
        double rp;          // radius of perigee of arc
        double a;           // semimajor axis of arc
        double e;           // eccentricity of arc
        double p;           // semilatus rectum arc
        double vv2 = norm(v2);
        ksi = 0.5 * vv2 * vv2 - (1.0 / norm(r1));
        h = norm(cross(r1, v1));
        p = h * h;
        a = -1.0 / (2 * a);
        e = sqrt((a - p) / a);
        rp = a * (1 - e);
        // Is the perigee lower than the radius of Earth?
        if (rp <= 1.0)
        {
            cout << "HIT!" << endl;
            return true;
        }
    }
    cout << "miss." << endl;
    return false;
}
