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
* $Id: Traj.cpp,v 1.16 2006/08/15 00:21:38 trs137 Exp $
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
        E(UNDEFINED),
        M(UNDEFINED),
        argLat(UNDEFINED),
        lonTrue(UNDEFINED),
        lonPer(UNDEFINED),
        r(0.0, 0.0, 0.0),
        v(0.0, 0.0, 0.0)
{
    // std::cout << "Traj constructor called with zero args.\n";
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
        E(UNDEFINED),
        M(UNDEFINED),
        argLat(UNDEFINED),
        lonTrue(UNDEFINED),
        lonPer(UNDEFINED),
        r(0.0, 0.0, 0.0),
        v(0.0, 0.0, 0.0)
{
    // std::cout << "Traj constructor called with classical elements.\n";
    randv(); // Calculate r and v vectors from the classical elements.
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
        E(UNDEFINED),
        M(UNDEFINED),
        argLat(UNDEFINED),
        lonTrue(UNDEFINED),
        lonPer(UNDEFINED),
        r(rin),
        v(vin)
{
    // std::cout << "Traj constructor called with state vector.\n";
    elorb();  // Calculate classical elements from the state vector.
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
        E(copy.E),
        M(copy.M),
        argLat(copy.argLat),
        lonTrue(copy.lonTrue),
        lonPer(copy.lonPer),
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
    E = t.E;
    M = t.M;
    argLat = t.argLat;
    lonTrue = t.lonTrue;
    lonPer = t.lonPer;
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
    cout << "TRAJECTORY PRINTOUT: " << endl;
    cout << "semimajor axis:               " << a << endl;
    cout << "eccentricity:                 " << e << endl;
    cout << "inclination:                  " << i << endl;
    cout << "RA of ascending node:         " << raan << endl;
    cout << "argument of periapsis:        " << w << endl;
    cout << "true anomaly:                 " << f << endl;
    cout << "eccentric/hyperbolic anomaly: " << E << endl;
    cout << "mean anomaly:                 " << M << endl;
    cout << "Argument of latitude          " << argLat << endl;
    cout << "True longitude                " << lonTrue << endl;
    cout << "Longitude of periapsis        " << lonPer << endl;
    cout << "radius vector:                " << r << ", magnitude = " << norm(r) << endl;
    cout << "velocity vector:              " << v << ", magnitude = " << norm(v) << endl;
}

/*
 *  Accessor methods for classic orbital elements.
 */
double Traj::get_a (void) { return a; }
double Traj::get_e (void) { return e; }
double Traj::get_i (void) { return i; }
double Traj::get_raan (void) { return raan; }
double Traj::get_w (void) { return w; }
double Traj::get_f (void) { return f; }
double Traj::get_E (void) { return E; }
double Traj::get_M (void) { return M; }
double Traj::get_argLat (void) { return argLat; }
double Traj::get_lonTrue (void) { return lonTrue; }
double Traj::get_lonPer (void) { return lonPer; }
Vec3 Traj::get_r (void) { return r; }
Vec3 Traj::get_v (void) { return v; }

/**
 * Mutator method for semimajor axis.
 * Constraints: [a < 0 : e > 1]
 *              [a = 0 : e = 1]
 *              [a > 0 : e < 1]
 */
void
Traj::set_a (double ain)
{
    if (    ((ain < 0) && (e <= 1))
         || ((fabs(ain) < SMALL) && (fabs(e-1.0) < SMALL))
         || ((ain > 0)  && (e >= 1)) )
    {
        cerr << "ERROR: Traj::set_a can't change semimajor axis to a value incompatible with eccentricity." << endl;
        exit(1);
    }
    a = ain;
}

/**
 * Mutator method for eccentricity.
 * Constraints: [a < 0 : e > 1]
 *              [a = 0 : e = 1]
 *              [a > 0 : e < 1]
 */
void
Traj::set_e (double ein)
{
    if (    ((a < 0) && (ein <= 1))
         || ((fabs(a) < SMALL) && (fabs(ein-1.0) < SMALL))
         || ((a > 0)  && (ein >= 1)) )
    {
        cerr << "ERROR: Traj::set_e can't change eccentricity to a value incompatible with semimajor axis." << endl;
        exit(1);
    }
    e = ein;
}
/**
 * Mutator method for inclination.
 * Constraints: [ 0 <= i <= 2*pi ]
 */
void
Traj::set_i (double iin)
{
    if ( (iin < 0) || (iin > 2*M_PI) )
    {
        cerr << "ERROR: Traj::set_i can't set bad inclination value." << endl;
        exit(1);
    }
    i = iin;
}

/**
 * Mutator method for RAAN.
 * Constraints: [ 0 <= RAAN <= 2*pi ]
 */
void
Traj::set_raan (double raanin)
{
    if ( (raanin < 0) || (raanin > 2*M_PI) )
    {
        cerr << "ERROR: Traj::set_raan can't set bad RAAN value." << endl;
        exit(1);
    }
    raan = raanin;
}

/**
 * Mutator method for argument of periapsis.
 * Constraints: [ 0 <= w <= 2*pi ]
 */
void
Traj::set_w (double win)
{
    if ( (win < 0) || (win > 2*M_PI) )
    {
        cerr << "ERROR: Traj::set_w can't set bad argument of periapsis value." << endl;
        exit(1);
    }
    w = win;
}

/**
 * Mutator method for true anomaly.
 * Constraints: [ 0 <= f <= 2*pi ]
 */
void
Traj::set_f (double fin)
{
    if ( (fin < 0) || (fin > 2*M_PI) )
    {
        cerr << "ERROR: Traj::set_f can't set bad true anomaly value." << endl;
        exit(1);
    }
    f = fin;
}

/**
 * Mutator method for geocentric position vector.
 * Constraints: [ norm(r) > 1.0 ] (units of earth radii)
 */
void
Traj::set_r (Vec3 rin)
{
    if ( norm(rin) < 1.0 )
    {
        cerr << "ERROR: Traj::set_r can't set position vector inside the Earth." << endl;
        exit(1);
    }
    r = rin;
}

/**
 * Mutator method for geocentric velocity vector.
 * Constraints: [ norm(v) < speed of light ] (units of ER/TU)
 */
void
Traj::set_v (Vec3 vin)
{
    if ( norm(vin) > 37922.655 )
    {
        cerr << "ERROR: I'm givin' it all she's got, Captain!" << endl;
        cerr << "ERROR: Traj::set_v can't set velocity vector faster than the speed of light." << endl;
        exit(1);
    }
    v = vin;
}

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
            cout << "NOTICE: circular equatorial orbit." << endl;
            f    = lonTrue; // set to True Longitude
        } else {
        // CIRCULAR INCLINED ORBIT
            w = 0.0;
            cout << "NOTICE: circular inclined orbit." << endl;
            f = argLat; // set to Argument of Latitude
        }
        else if (( i < SMALL) || (fabs(i-M_PI)) < SMALL)
        {
        // ELLIPTICAL EQUATORIAL
            cout << "NOTICE: noncircular equatorial traj." << endl;
            w = lonPer; // set to Longitude of Periapsis
            raan = 0.0;
        }

    // Calculate semilatus rectum or semiparameter
    double p = a * (1 - e*e);

    // Calculate some temporary values
    double cos_f = cos(f);
    double sin_f = sin(f);
    double temp  = p / (1.0 + e * cos_f);

    // Calculate position and velocity in PQW frame.
    Vec3 r_pqw((temp*cos_f), (temp*sin_f), (0.0));
    if (fabs(p) < SMALL) p = SMALL;
    Vec3 v_pqw((-sin_f/sqrt(p)), ((e+cos_f)/sqrt(p)), 0.0);

    // Transform PQW to Geocentric Equitorial
    // r = r_pqw;
    // v = v_pqw;

    r = rotZ(rotX(rotZ(r_pqw, w), i), raan);
    v = rotZ(rotX(rotZ(v_pqw, w), i), raan);

    // Now calculte the miscellaneous stuff.
    double vv = norm(v);
    double rr = norm(r);
    Vec3 h = cross(r, v);
    Vec3 nodeVector = cross(Vec3(0,0,1), h);  // nodeVector = the node vector
    double nn = norm(nodeVector);

    // Argument of Latitude (for circular inclined orbits)
    if (e < SMALL)
    {
        argLat = acos(dot(nodeVector, r) / (nn * rr));
        // Quadrant check!
        if (r.getZ() < 0) argLat = 2*M_PI - argLat;
    }else{
        argLat = UNDEFINED;
    }

    // True Longitude (for circular equatorial orbits)
    if ((e < SMALL) && (i < SMALL))
    {
        lonTrue = acos(r.getX() / rr );
        // Quadrant check!
        if (r.getY() < 0 ) lonTrue = 2*M_PI - lonTrue;
    } else {
        lonTrue = UNDEFINED;
    }

    // Longitude of Periapsis (for non-circular equatorial trajs)
    Vec3 eccVector = ((vv*vv - 1.0/rr)*r - dot(r,v)*v); // CANONICAL UNITS ONLY
    if (i < SMALL)
    {
        lonPer = acos( eccVector.getX() / e );
        if (eccVector.getY() < 0 ) lonPer = 2*M_PI - lonPer;
    } else {
        lonPer = UNDEFINED;
    }

    // Eccentric Anomaly
    double num = tan(f/2);
    double denom = sqrt( (1+e)/(1-e) );
    E = 2 * atan2 (num, denom);

    // Mean Anomaly
    M = E - e * sin(E);
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

    // Argument of Latitude (for circular inclined orbits)
    if (e < SMALL)
    {
        argLat = acos(dot(nodeVector, r) / (nn * rr));
        // Quadrant check!
        if (r.getZ() < 0) argLat = 2*M_PI - argLat;
    }else{
        argLat = UNDEFINED;
    }

    // True Longitude (for circular equatorial orbits)
    if ((e < SMALL) && (i < SMALL))
    {
        lonTrue = acos(r.getX() / rr );
        // Quadrant check!
        if (r.getY() < 0 ) lonTrue = 2*M_PI - lonTrue;
    } else {
        lonTrue = UNDEFINED;
    }

    // Longitude of Periapsis (for non-circular equatorial trajs)
    if (i < SMALL)
    {
        lonPer = acos( eccVector.getX() / e );
        if (eccVector.getY() < 0 ) lonPer = 2*M_PI - lonPer;
    } else {
        lonPer = UNDEFINED;
    }

    // Eccentric Anomaly
    double num = tan(f/2);
    double denom = sqrt( (1+e)/(1-e) );
    E = 2 * atan2 (num, denom);

    // Mean Anomaly
    M = E - e * sin(E);
}
