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
* $Id: Elements.h,v 1.2 2006/08/06 00:33:48 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

/** @file
 * Common functions related to classical orbital elements. 
 */

#ifndef _ELEMENTS_H_
#define _ELEMENTS_H_

#include <math.h>
#include "Orbgnosis.h"
#include "Traj.h"
#include "Vec3.h"

/**
 * Returns classical orbital elements given radius and velocity vectors
 * in canonical units.
 */
Traj
elorb(Vec3 r, Vec3 v)
{
    double vv = norm(v);
    double rr = norm(r);
    Vec3 h = cross(r, v);
    double hh = norm(h);   // h = specific angular momentum.  Units ER^2 / TU
    Vec3 nodeVector = cross(Vec3(0,0,1), h);  // nodeVector = the node vector
    double nn = norm(nodeVector);

    // Eccentricity
    Vec3 e = ((vv*vv - 1.0/rr)*r - dot(r,v)*v); // CANONICAL UNITS ONLY
    // Otherwise (S.I. units, MU = 398600 km^3 / s^2)
    // e = (1/MU) * ((vv*vv - MU/rr)*r - dot(r,v)*v);
    double ee = norm(e);

    // Specific Mechanical Energy, greek letter ksi
    // Units: ER^2 / TU^2
    double ksi = vv*vv/2.0 - 1.0/rr; // CANONICAL UNITS ONLY
    // double ksi = vv*vv/2.0 - MU/rr // general for non-canonical

    // Semimajor axis.  Units: ER
    double a = -1.0 / (2.0 * ksi);
    // double a = -MU / (2.0 * ksi); // general for non-canonical

    // Inclination.  Units: Radians
    // Quadrant check is not necessary.
    double i = acos(h.getZ() / hh);

    // RA of ascending node.  Units: Radians.
    double raan = acos(nodeVector.getX() / nn);
    // Quadrant check!
    if (nodeVector.getY() < 0 ) raan = 2*M_PI - raan;

    // Argument of Perigee. Units: Radians.
    double w = acos( dot(nodeVector, e) / (nn * ee));
    // Quadrant check!
    if (e.getZ() < 0) w = 2*M_PI - w;

    // True Anomaly.  Units: Radians.
    double f = acos( dot(e, r) / (ee * rr));
    // Quadrant check!
    if (dot(r,v) < 0) f = 2*M_PI - f;

    return Traj(a, ee, i, raan, w, f);
}

#endif /* _ELEMENTS_H_ */
