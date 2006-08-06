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
* $Id: Traj.cpp,v 1.12 2006/08/06 22:36:13 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

#include <math.h>
#include "Orbgnosis.h"
#include "Traj.h"
#include <iostream>
using namespace std;

/**
 * The Traj constructor.
 */
Traj::Traj (void) :
        a(0.0),
        e(0.0),
        i(0.0),
        raan(0.0),
        w(0.0),
        f(0.0),
        r(0.0, 0.0, 0.0),
        v(0.0, 0.0, 0.0)
{
    std::cout << "Traj constructor called with zero args.\n";
}

/**
 * Traj constructor.
 * All six classical Keplerian elements are
 * specified in the constructor.
 * @param ain the semimajor axis (length units).
 * @param ein the eccentricity (dimensionless).
 * @param iin the inclination (radians).
 * @param raanin the right ascention of the ascending node (radians)
 * @param win the argument of periapsis (radians).
 * @param fin the true anomaly (radians).
 */
Traj::Traj (double ain, double ein, double iin, double raanin,
            double win, double fin) :
        a(ain),
        e(ein),
        i(iin),
        raan(raanin),
        w(win),
        f(fin),
        r(0.0, 0.0, 0.0),
        v(0.0, 0.0, 0.0)
{
    std::cout << "Traj constructor called with classical elements.\n";
    randv();
}

/**
 * Traj constructor.
 * State vector in canonical units used.
 */
Traj::Traj (Vec3 rin, Vec3 vin) :
        a(0.0),
        e(0.0),
        i(0.0),
        raan(0.0),
        w(0.0),
        f(0.0),
        r(rin),
        v(vin)
{
    std::cout << "Traj constructor called with state vector.\n";
    elorb();
}

/**
 * The Traj copy constructor.
 * @param copy the reference to the Traj to be copied.
 */
Traj::Traj (const Traj& copy) :
        a(copy.a),
        e(copy.e),
        i(copy.i),
        raan(copy.raan),
        w(copy.w),
        f(copy.f),
        r(copy.r),
        v(copy.v)
{
    //std::cout << "Traj copy constructor called.\n";
}

/**
 * The Traj copy assignment operator.
 * @param t the Traj on the right hand side of the equality.
 */
Traj&
Traj::operator = (Traj t)
{
    a = t.a;
    e = t.e;
    i = t.i;
    raan = t.raan;
    w = t.w;
    f = t.f;
    r = t.r;
    v = t.v;
    return *this;
}

/**
 * The Traj destructor.
 */
Traj::~Traj (void)
{
    // cout << "Traj destructor called.\n";
}

/**
 * Prints classic orbital elements nicely.
 */
void
Traj::print (void)
{
    cout << "semimajor axis:        " << a << endl;
    cout << "eccentricity:          " << e << endl;
    cout << "inclination:           " << i << endl;
    cout << "RA of ascending node:  " << raan << endl;
    cout << "argument of periapsis: " << w << endl;
    cout << "true anomaly:          " << f << endl;
    cout << "radius vector:         " << r << ", " << norm(r) << endl;
    cout << "velocity vector:       " << v << ", " << norm(v) << endl;
}

/*
 *  Accessor & Mutator methods for classic orbital elements.
 */
double Traj::get_a (void) { return a; }
double Traj::get_e (void) { return e; }
double Traj::get_i (void) { return i; }
double Traj::get_raan (void) { return raan; }
double Traj::get_w (void) { return w; }
double Traj::get_f (void) { return f; }
Vec3 Traj::get_r (void) { return r; }
Vec3 Traj::get_v (void) { return v; }
void Traj::set_a (double ain) { a = ain; }
void Traj::set_e (double ein) { a = ein; }
void Traj::set_i (double iin) { a = iin; }
void Traj::set_raan (double raanin) { a = raanin; }
void Traj::set_w (double win) { a = win; }
void Traj::set_f (double fin) { a = fin; }
void Traj::set_r (Vec3 rin) { r = rin; }
void Traj::set_v (Vec3 vin) { v = vin; }

/**
 * Calculates position and velocity vectors, given
 * classical orbital elements in canonical units.
 */
void
Traj::randv()
{
    // Set up angles for certain special case orbits.
    //
    if (e < SMALL)
        if ((i < SMALL) || (fabs(i-M_PI)) < SMALL)
        // CIRCULAR EQUATORIAL ORBIT
        {
            w    = 0.0; 
            raan = 0.0;
            cout << "ERROR: circular equatorial orbit." << endl;
            // f    = FIXME;  // set to True Longitude
        } else {
        // CIRCULAR INCLINED ORBIT
            w = 0.0;
            cout << "ERROR: circular inclined orbit." << endl;
            // f = FIXME; // set to Argument of Latitude
        }
        else if (( i < SMALL) || (fabs(i-M_PI)) < SMALL)
        {
        // ELLIPTICAL EQUATORIAL
            cout << "ERROR: elliptical equatorial orbit." << endl;
            // w = FIXME; // set to Longitude of Periapsis
            raan = 0.0;
        }

    // Calculate semilatus rectum or semiparameter
    double p = a * (1 - e*e);

    // Calculate some temporary values
    double cos_f = cos(f);
    double sin_f = sin(f);
    double temp  = p / (1.0 + e * cos_f);

    Vec3 rPQW((temp*cos_f), (temp*sin_f), (0.0));
    if (fabs(p) < SMALL) p = SMALL;
    Vec3 vPQW((-sin_f/sqrt(p)), ((e+cos_f)/sqrt(p)), 0.0 );

    // Transform PQW to Geocentric Equitorial

    Vec3 TempVec;

    TempVec = rotZ(rPQW, -w);
    TempVec = rotX(TempVec, -i);
    r = rotZ(TempVec, -raan);

    TempVec = rotZ(vPQW, -w);
    TempVec = rotX(TempVec, -i);
    v = rotZ(TempVec, -raan);
}




/**
 * Calculates classical orbital elements, given radius and velocity
 * in canonical units.
 */
void
Traj::elorb(void)
{
    double vv = norm(v);
    double rr = norm(r);
    Vec3 h = cross(r, v);
    double hh = norm(h);   // h = specific angular momentum.  Units ER^2 / TU
    Vec3 nodeVector = cross(Vec3(0,0,1), h);  // nodeVector = the node vector
    double nn = norm(nodeVector);

    // Specific Mechanical Energy, greek letter ksi.  Units: ER^2 / TU^2
    double ksi = vv*vv/2.0 - 1.0/rr; // CANONICAL UNITS ONLY
    // double ksi = vv*vv/2.0 - MU/rr // general for non-canonical

    // Semimajor axis.  Units: ER
    a = -1.0 / (2.0 * ksi);
    // a = -MU / (2.0 * ksi); // general for non-canonical

    // Eccentricity.  Unitless.
    Vec3 eccVector = ((vv*vv - 1.0/rr)*r - dot(r,v)*v); // CANONICAL UNITS ONLY
    // Otherwise (S.I. units, MU = 398600 km^3 / s^2)
    // eVector = (1/MU) * ((vv*vv - MU/rr)*r - dot(r,v)*v);
    e = norm(eccVector);

    // Inclination.  Units: Radians
    // Quadrant check is not necessary.
    i = acos(h.getZ() / hh);

    // RA of ascending node.  Units: Radians.
    raan = acos(nodeVector.getX() / nn);
    // Quadrant check!
    if (nodeVector.getY() < 0 ) raan = 2*M_PI - raan;

    // Argument of Perigee. Units: Radians.
    w = acos( dot(nodeVector, eccVector) / (nn * e));
    // Quadrant check!
    if (eccVector.getZ() < 0) w = 2*M_PI - w;

    // True Anomaly.  Units: Radians.
    f = acos( dot(eccVector, r) / (e * rr));
    // Quadrant check!
    if (dot(r,v) < 0) f = 2*M_PI - f;
}

