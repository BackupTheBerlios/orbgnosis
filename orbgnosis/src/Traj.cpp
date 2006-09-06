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
* $Id: Traj.cpp,v 1.20 2006/09/06 17:34:26 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

#include "Orbgnosis.h"
#include "Traj.h"
#include "Vec3.h"
#include <iostream>
#include <math.h>
using namespace std;

/**
 * The Traj constructor.
 */
Traj::Traj (void) : a(0.0), e(0.0), i(0.0), raan(0.0), w(0.0), f(0.0),
    r(0.0, 0.0, 0.0), v(0.0, 0.0, 0.0),
    E(NAN), M(NAN), argLat(NAN), lonTrue(NAN), lonPer(NAN),
    e_vector(0.0, 0.0, 0.0),
    h_vector(0.0, 0.0, 0.0),
    n_vector(0.0, 0.0, 0.0)
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
Traj::Traj (double ain, double ein, double iin, double raanin, double win,
    double fin) :
    a(ain), e(ein), i(iin), raan(raanin), w(win), f(fin),
    r(0.0, 0.0, 0.0), v(0.0, 0.0, 0.0),
    E(NAN), M(NAN), argLat(NAN), lonTrue(NAN), lonPer(NAN),
    e_vector(0.0, 0.0, 0.0),
    h_vector(0.0, 0.0, 0.0),
    n_vector(0.0, 0.0, 0.0)
{
    // std::cout << "Traj constructor called with classical elements.\n";
    randv(); // Calculate r and v vectors from the classical elements.
}

/**
 * Traj constructor.
 * State vector in canonical units used.
 */
Traj::Traj (Vec3 rin, Vec3 vin) : 
    a(0.0), e(0.0), i(0.0), raan(0.0), w(0.0), f(0.0),
    r(rin), v(vin),
    E(NAN), M(NAN), argLat(NAN), lonTrue(NAN), lonPer(NAN),
    e_vector(0.0, 0.0, 0.0),
    h_vector(0.0, 0.0, 0.0),
    n_vector(0.0, 0.0, 0.0)
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
        r(copy.r),
        v(copy.v),
        E(copy.E),
        M(copy.M),
        argLat(copy.argLat),
        lonTrue(copy.lonTrue),
        lonPer(copy.lonPer),
        e_vector(copy.e_vector),
        h_vector(copy.h_vector),
        n_vector(copy.n_vector)
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
    if (this != &t) // avoid pointless self-assignment
    {
        a = t.a;
        e = t.e;
        i = t.i;
        raan = t.raan;
        w = t.w;
        f = t.f;
        r = t.r;
        v = t.v;
        E = t.E;
        M = t.M;
        argLat = t.argLat;
        lonTrue = t.lonTrue;
        lonPer = t.lonPer;
        e_vector = t.e_vector;
        h_vector = t.h_vector;
        n_vector = t.n_vector;
    }
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
    cout << "eccentricity:                 " << e_vector << ", " << e << endl;
    cout << "inclination:                  " << i << endl;
    cout << "RA of ascending node:         " << raan << endl;
    cout << "argument of periapsis:        " << w << endl;
    cout << "true anomaly:                 " << f << endl;
    cout << "radius vector:                " << r << ", norm = " << norm(r) << endl;
    cout << "velocity vector:              " << v << ", norm = " << norm(v) << endl;
    cout << "eccentric/hyperbolic anomaly: " << E << endl;
    cout << "mean anomaly:                 " << M << endl;
    cout << "Argument of latitude          " << argLat << endl;
    cout << "True longitude                " << lonTrue << endl;
    cout << "Longitude of periapsis        " << lonPer << endl;
    cout << "spec angular momentum:        " << h_vector << ", norm = " << norm(h_vector) << endl;
    cout << "node vector:                  " << n_vector << ", norm = " << norm(n_vector) << endl;
}

/**
 * Solve Kepler's problem.  Given a state vector (Traj) and a time interval,
 * find the state vector after the time interval has elapsed.
 */
Traj
kepler (Traj traj_0, double t)
{
    if (fabs(t) <= SMALL) return traj_0;
    else {
    // set up local variables
    double r0;      // initial radius
    double v0;      // initial velocity
    double F, G;    // universal variable f and g expressions
    double Fdot, Gdot; // F and G's rates of change
    double Xold;    // universal variable
    double Xnew;    // universal variable
    double Xold2;   // Xold squared
    double Xnew2;   // Xnew squared
    double Znew;    // new value of Z
    double C2new;   // Stumpff C2
    double C3new;   // Stumpff C3
    double tnew;    // new time
    double rdotv;   // r0 dot v0;
    double a;       // semimajor axis;
    double alpha;   // 1 / (semimajor axis)
    double ksi;     // specific mechanicak energy
    double period;  // orbital period
    double S, W;    // variables for parabolic special case
        
    r0 = norm(traj_0.get_r()); // current radius
    v0 = norm(traj_0.get_v()); // current velocity
    Xold = 0.0;
    Znew = 0.0;
    rdotv = dot(traj_0.get_r(), traj_0.get_v());

    ksi = (0.5 * v0 * v0) - (1 / r0);     // canonical!
    alpha = -2.0 * ksi;

    a = traj_0.get_a();

    if (alpha >= SMALL)
    {
        period = 2*M_PI * sqrt(pow(fabs(a), 3.0));
        if (fabs(t) > fabs(period)) t = fmod (t, period); // multirev
        if (fabs(alpha - 1.0) > SMALL) Xold = t * alpha;
        else
             Xold = t * alpha * 0.97;
    } else {
        if (fabs(alpha) < SMALL)
        {
            // Parabola XXX TESTME
            double h = norm(traj_0.get_h_vector());
            double p = h * h;
            S = 0.5 * (M_PI/2.0 - atan(3.0 * sqrt(1.0 / (p * p * p)) * t));
            W = atan(pow(tan(S), 1.0 / 3.0));
            Xold = sqrt(p) * (2.0 * (1.0/tan(2.0 * W)) );
            alpha = 0.0;
        } else {
            // Hyperbola TODO
            // This only works correctly for positive t.
            double temp = -2.0 * t / 
                (a * (rdotv + sqrt(-a) * (1.0 - r0 * alpha)));
            Xold = sqrt(-a) * log(temp);
        }
    }

    return traj_t; }
}

/*
 *  Accessors.
 */
double Traj::get_a (void) { return a; }
double Traj::get_e (void) { return e; }
double Traj::get_i (void) { return i; }
double Traj::get_raan (void) { return raan; }
double Traj::get_w (void) { return w; }
double Traj::get_f (void) { return f; }
Vec3 Traj::get_r (void) { return r; }
Vec3 Traj::get_v (void) { return v; }
double Traj::get_E (void) { return E; }
double Traj::get_M (void) { return M; }
double Traj::get_argLat (void) { return argLat; }
double Traj::get_lonTrue (void) { return lonTrue; }
double Traj::get_lonPer (void) { return lonPer; }
Vec3 Traj::get_e_vector (void) { return e_vector; }
Vec3 Traj::get_h_vector (void) { return h_vector; }
Vec3 Traj::get_n_vector (void) { return n_vector; }

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
    randv(); // changing classical orbital element requires re-running randv()
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
    randv(); // changing classical orbital element requires re-running randv()
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
    randv(); // changing classical orbital element requires re-running randv()
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
    randv(); // changing classical orbital element requires re-running randv()
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
    randv(); // changing classical orbital element requires re-running randv()
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
    randv(); // changing classical orbital element requires re-running randv()
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
    elorb();  // changing state vector requires re-running elorb()
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
    elorb();  // changing state vector requires re-running elorb()
}

/**
 * Calculates position and velocity vectors, given
 * classical orbital elements in canonical units.
 */
void
Traj::randv()
{
    // Set up angles for certain special case orbits.
    // This assumes that if special non-classical orbital
    // Elements are needed to specify the orbit they will
    // have already been defined.

    // Calculate semilatus rectum or semiparameter
    double p = a * (1 - e*e);

    // Calculate some temporary values
    double cos_f = cos(f);
    double sin_f = sin(f);
    double temp  = p / (1.0 + e * cos_f);

    // Calculate position and velocity in PQW frame of reference.
    Vec3 r_pqw((temp*cos_f), (temp*sin_f), (0.0));
    if (fabs(p) < SMALL) p = SMALL;
    Vec3 v_pqw((-sin_f/sqrt(p)), ((e+cos_f)/sqrt(p)), 0.0);

    // Transform orbital frame of reference to geocentric equitorial
    r = rotZ(rotX(rotZ(r_pqw, w), i), raan);
    v = rotZ(rotX(rotZ(v_pqw, w), i), raan);

    // Fill in e_ h_ and n_vectors.
    double rr = norm(r);
    double vv = norm(v);
    e_vector = ((vv*vv - 1.0/rr)*r - dot(r,v)*v); // CANONICAL UNITS ONLY
    h_vector = cross(r, v);
    n_vector = cross(Vec3(0,0,1), h_vector);

    // Everything is solved except for the other anomalies
    // and special-case orbital elements.
    // randv and elorb share this stuff, so it has its own function.
    anomalies();
    special();
} // end randv

/**
 * Calculates classical orbital elements, given radius and velocity
 * in canonical units.
 */
void
Traj::elorb(void)
{
    h_vector = cross(r, v);
    n_vector = cross(Vec3(0,0,1), h_vector);
    double vv = norm(v);
    double rr = norm(r);
    double hh = norm(h_vector);
    double nn = norm(n_vector);

    // Specific Mechanical Energy, greek letter ksi.  Units: ER^2 / TU^2
    double ksi = vv*vv/2.0 - 1.0/rr; // CANONICAL UNITS ONLY
    // double ksi = vv*vv/2.0 - MU/rr // general for non-canonical

    // Semimajor axis.  Units: ER
    a = -1.0 / (2.0 * ksi);
    // a = -MU / (2.0 * ksi); // general for non-canonical

    // Eccentricity.  Unitless.
    e_vector = ((vv*vv - 1.0/rr)*r - dot(r,v)*v); // CANONICAL UNITS ONLY
    // Otherwise (S.I. units, MU = 398600 km^3 / s^2)
    // eVector = (1/MU) * ((vv*vv - MU/rr)*r - dot(r,v)*v);
    e = norm(e_vector);

    // Inclination.  Units: Radians
    // Quadrant check is not necessary.
    i = acos(h_vector.getZ() / hh);

    // RA of ascending node.  Units: Radians.
    raan = acos(n_vector.getX() / nn);
    // Quadrant check!
    if (n_vector.getY() < 0 ) raan = 2*M_PI - raan;

    // Argument of Perigee. Units: Radians.
    w = acos( dot(n_vector, e_vector) / (nn * e));
    // Quadrant check!
    if (e_vector.getZ() < 0) w = 2*M_PI - w;

    // True Anomaly.  Units: Radians.
    f = acos( dot(e_vector, r) / (e * rr));
    // Quadrant check!
    if (dot(r,v) < 0) f = 2*M_PI - f;

    // Everything is solved except for the other anomalies
    // and special-case orbital elements.
    anomalies();
    special();
}

void
Traj::anomalies()
{
    if ( e < 1.0 ) // elliptical
    {
        double num = tan(f/2);
        double denom = sqrt( (1+e)/(1-e) );
        E = 2 * atan2 (num, denom);
        if (E < 0.0) E = 2*M_PI + E;
        M = E - e * sin(E);
    } else if (fabs(e - 1.0) < SMALL) { // parabolic
        E = tan(f/2); // hyperbolic anomaly is simple
    } else { // hyperbolic 
        E = 2 * atanh ( sqrt( (e-1)/(e+1) ) * tan(f/2)) ;
        M = e * sinh(E) - E;
    }
    argLat = w + f;
}

void
Traj::special()
{
    double rr = norm(r);
    // Argument of Latitude (only for circular inclined orbits)
    if (e < SMALL)
    {
        argLat = acos(dot(n_vector, r) / (norm(n_vector) * rr));
        // Quadrant check!
        if (r.getZ() < 0) argLat = 2*M_PI - argLat;
    }

    // True Longitude (only for circular equatorial orbits)
    if ((e < SMALL) && (i < SMALL))
    {
        lonTrue = acos(r.getX() / rr );
        // Quadrant check!
        if (r.getY() < 0 ) lonTrue = 2*M_PI - lonTrue;
    }

    // Longitude of Periapsis (only for non-circular equatorial trajs)
    if (i < SMALL)
    {
        lonPer = acos( e_vector.getX() / e );
        if (e_vector.getY() < 0 ) lonPer = 2*M_PI - lonPer;
    }
}
