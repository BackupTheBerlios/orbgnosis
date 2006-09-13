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
* $Id: Traj.cpp,v 1.26 2006/09/13 02:01:15 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

#include "Orbgnosis.h"
#include "Stumpff.h"
#include "Traj.h"
#include "Vec3.h"
#include <iostream>
#include <math.h>

using namespace std;

/**
 * The Traj constructor.
 */
Traj::Traj ( void ) : a( 0.0 ), e( 0.0 ), i( 0.0 ), raan( 0.0 ), w( 0.0 ), f( 0.0 ),
        r( 0.0, 0.0, 0.0 ), v( 0.0, 0.0, 0.0 ),
        E( NAN ), M( NAN ), argLat( NAN ), lonTrue( NAN ), lonPer( NAN ),
        e_vector( 0.0, 0.0, 0.0 ),
        h_vector( 0.0, 0.0, 0.0 ),
        n_vector( 0.0, 0.0, 0.0 )
{
    // std::cout << "Traj constructor called with zero args.\n";
}

/**
 * Traj constructor with classical elements.
 * The six classical Keplerian elements are specified in this version
 * of the constructor.  The state vectors (position and velocity) are
 * calculated automatically.
 * @param ain semimajor axis (length units).
 * @param ein eccentricity (dimensionless).
 * @param inclination (radians).
 * @param raanin right ascention of the ascending node (radians)
 * @param win argument of periapsis (radians).
 * @param fin true anomaly (radians).
 */
Traj::Traj ( double ain, double ein, double iin, double raanin, double win,
             double fin ) :
        a( ain ), e( ein ), i( iin ), raan( raanin ), w( win ), f( fin ),
        r( 0.0, 0.0, 0.0 ), v( 0.0, 0.0, 0.0 ),
        E( NAN ), M( NAN ), argLat( NAN ), lonTrue( NAN ), lonPer( NAN ),
        e_vector( 0.0, 0.0, 0.0 ),
        h_vector( 0.0, 0.0, 0.0 ),
        n_vector( 0.0, 0.0, 0.0 )
{
    // std::cout << "Traj constructor called with classical elements.\n";
    randv(); // Calculate r and v vectors from the classical elements.
}

/**
 * Traj constructor with state vector.
 * This version of the constructor takes position and velocity vectors
 * in canonical units.  Classical elements are calculated automatically from
 * these vectors.
 * @param rin radius in geocentric IJK frame (Earth radii)
 * @param vin velocity in geocentric IJK frame (Earth radii per canonical time unit)
 */
Traj::Traj ( Vec3 rin, Vec3 vin ) :
        a( 0.0 ), e( 0.0 ), i( 0.0 ), raan( 0.0 ), w( 0.0 ), f( 0.0 ),
        r( rin ), v( vin ),
        E( NAN ), M( NAN ), argLat( NAN ), lonTrue( NAN ), lonPer( NAN ),
        e_vector( 0.0, 0.0, 0.0 ),
        h_vector( 0.0, 0.0, 0.0 ),
        n_vector( 0.0, 0.0, 0.0 )
{
    // std::cout << "Traj constructor called with state vector.\n";
    elorb();  // Calculate classical elements from the state vector.
}

/**
 * The Traj copy constructor.
 * @param copy the reference to the Traj to be copied.
 */
Traj::Traj ( const Traj& copy ) :
        a( copy.a ),
        e( copy.e ),
        i( copy.i ),
        raan( copy.raan ),
        w( copy.w ),
        f( copy.f ),
        r( copy.r ),
        v( copy.v ),
        E( copy.E ),
        M( copy.M ),
        argLat( copy.argLat ),
        lonTrue( copy.lonTrue ),
        lonPer( copy.lonPer ),
        e_vector( copy.e_vector ),
        h_vector( copy.h_vector ),
        n_vector( copy.n_vector )
{
    //std::cout << "Traj copy constructor called.\n";
}

/**
 * The Traj copy assignment operator.
 * @param t the Traj on the right hand side of the equality.
 */
Traj&
Traj::operator = ( Traj t )
{
    if ( this != &t )   // avoid pointless self-assignment
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
Traj::~Traj ( void )
{
    // cout << "Traj destructor called.\n";
}

/**
 * One way to print a trajectory's data.
 */
void
Traj::print ( void )
{
    cout << "TRAJECTORY PRINTOUT: " << endl;
    cout << "semimajor axis:               " << a << endl;
    cout << "eccentricity:                 " << e << ", " << e_vector << endl;
    cout << "inclination:                  " << i << endl;
    cout << "RA of ascending node:         " << raan << endl;
    cout << "argument of periapsis:        " << w << endl;
    cout << "true anomaly:                 " << f << endl;
    cout << "radius vector:                " << r << ", norm = " << norm( r ) << endl;
    cout << "velocity vector:              " << v << ", norm = " << norm( v ) << endl;
    cout << "eccentric/hyperbolic anomaly: " << E << endl;
    cout << "mean anomaly:                 " << M << endl;
    cout << "Argument of latitude          " << argLat << endl;
    cout << "True longitude                " << lonTrue << endl;
    cout << "Longitude of periapsis        " << lonPer << endl;
    cout << "spec angular momentum:        " << h_vector << ", norm = " << norm( h_vector ) << endl;
    cout << "node vector:                  " << n_vector << ", norm = " << norm( n_vector ) << endl;
}

/*
 *  Accessors.
 */
double Traj::get_a ( void )
{
    return a;
}

double Traj::get_e ( void )
{
    return e;
}

double Traj::get_i ( void )
{
    return i;
}

double Traj::get_raan ( void )
{
    return raan;
}

double Traj::get_w ( void )
{
    return w;
}

double Traj::get_f ( void )
{
    return f;
}

Vec3 Traj::get_r ( void )
{
    return r;
}

Vec3 Traj::get_v ( void )
{
    return v;
}

double Traj::get_E ( void )
{
    return E;
}

double Traj::get_M ( void )
{
    return M;
}

double Traj::get_argLat ( void )
{
    return argLat;
}

double Traj::get_lonTrue ( void )
{
    return lonTrue;
}

double Traj::get_lonPer ( void )
{
    return lonPer;
}

Vec3 Traj::get_e_vector ( void )
{
    return e_vector;
}

Vec3 Traj::get_h_vector ( void )
{
    return h_vector;
}

Vec3 Traj::get_n_vector ( void )
{
    return n_vector;
}

/**
 * Completely re-set a trajectory by specifying all six classical orbital
 * elements.  This method automatically takes care of the state vector
 * and miscelleous things.
 */
void
Traj::set_elorb ( double ain, double ein, double iin, double raanin, double win,
                  double fin)
{
    a = ain;
    e = ein;
    i = iin;
    raan = raanin;
    w = win;
    f = fin;
    randv();
}

/**
 * Completely re-set a traj by specifying the position and velocity vectors.
 * This method automatically takes care of the classical elements and misc things.
 */
void
Traj::set_randv ( Vec3 rin, Vec3 vin )
{
    r = rin;
    v = vin;
    elorb();
}

/**
 * Mutator method for semimajor axis.  Automatically re-calculates state vector.
 * @param ain the semimajor axis must correspond appropriately to the
 * eccentricity.  For e > 1, ain must be negative.  For e = 1, ain must be zero.
 * For e < 1, ain must be positive.
 */
void
Traj::set_a ( double ain )
{
    if ( ( ( ain < 0 ) && ( e <= 1 ) )
            || ( ( fabs( ain ) < SMALL ) && ( fabs( e - 1.0 ) < SMALL ) )
            || ( ( ain > 0 ) && ( e >= 1 ) ) )
    {
        cerr << "ERROR: Traj::set_a can't change semimajor axis to a value incompatible with eccentricity." << endl;
        exit( 1 );
    }

    a = ain;
    randv(); // changing classical orbital element requires re-running randv()
}

/**
 * Mutator method for eccentricity.  Automatically re-calculates state vector.
 * @param ein the eccentricity must correspond correctly to the
 * semimajor axis, a.  For negative a, ein must be greater than one.
 * For a = 0, e must be exactly one.  For positive a, e must be between
 * zero and one, non-inclusive.
 */
void
Traj::set_e ( double ein )
{
    if ( ( ( a < 0 ) && ( ein <= 1 ) )
            || ( ( fabs( a ) < SMALL ) && ( fabs( ein - 1.0 ) < SMALL ) )
            || ( ( a > 0 ) && ( ( ein >= 1 ) || ( ein <= 0 ) ) ) )
    {
        cerr << "ERROR: Traj::set_e can't change eccentricity to a value incompatible with semimajor axis." << endl;
        exit( 1 );
    }

    e = ein;
    randv(); // changing classical orbital element requires re-running randv()
}

/**
 * Mutator method for inclination.  Automatically re-calculates state vector.
 * @param iin the inclination must be between zero and two pi, inclusive.
 */
void
Traj::set_i ( double iin )
{
    if ( ( iin < 0 ) || ( iin > 2 * M_PI ) )
    {
        cerr << "ERROR: Traj::set_i can't set bad inclination value." << endl;
        exit( 1 );
    }

    i = iin;
    randv(); // changing classical orbital element requires re-running randv()
}

/**
 * Mutator method for RAAN.  Automatically re-calculates state vector.
 * @param raanin the right ascention of the ascending node must be between zero
 * and two pi, inclusive.
 * Constraints: [ 0 <= RAAN <= 2*pi ]
 */
void
Traj::set_raan ( double raanin )
{
    if ( ( raanin < 0 ) || ( raanin > 2 * M_PI ) )
    {
        cerr << "ERROR: Traj::set_raan can't set bad RAAN value." << endl;
        exit( 1 );
    }

    raan = raanin;
    randv(); // changing classical orbital element requires re-running randv()
}

/**
 * Mutator method for argument of periapsis. Automatically re-calculates state vector.
 * @param win the argument of periapsis must bet between zero and two pi, inclusive.
 */
void
Traj::set_w ( double win )
{
    if ( ( win < 0 ) || ( win > 2 * M_PI ) )
    {
        cerr << "ERROR: Traj::set_w can't set bad argument of periapsis value." << endl;
        exit( 1 );
    }

    w = win;
    randv(); // changing classical orbital element requires re-running randv()
}

/**
 * Mutator method for true anomaly. Automatically re-calculates state vector.
 * @param f the true anomaly must be between zero and two pi, inclusive.
 */
void
Traj::set_f ( double fin )
{
    if ( ( fin < 0 ) || ( fin > 2 * M_PI ) )
    {
        cerr << "ERROR: Traj::set_f can't set to bad true anomaly value." << endl;
        exit( 1 );
    }

    f = fin;
    randv(); // changing classical orbital element requires re-running randv()
}

/**
 * Mutator method for geocentric position vector.  Automatically re-calculates
 * classical elements.
 * @param rin the geocentric position vector must not be inside the Earth.
 * Units of length are in Earth radii.
 */
void
Traj::set_r ( Vec3 rin )
{
    if ( norm( rin ) < 1.0 )
    {
        cerr << "ERROR: Traj::set_r can't set position vector inside the Earth." << endl;
        exit( 1 );
    }

    r = rin;
    elorb();  // changing state vector requires re-running elorb()
}

/**
 * Mutator method for geocentric velocity vector. Automatically re-calculates
 * classical elements.
 * @param vin the geocentric velocity vector must not be faster than the speed
 * of light.  Units are ER/TU.
 */
void
Traj::set_v ( Vec3 vin )
{
    if ( norm( vin ) > 37922.655 )
    {
        cerr << "ERROR: I'm givin' it all she's got, Captain!" << endl;
        cerr << "ERROR: Traj::set_v can't set velocity vector faster than the speed of light." << endl;
        exit( 1 );
    }

    v = vin;
    elorb();  // changing state vector requires re-running elorb()
}

/**
 * Calculates position and velocity vectors, given
 * classical orbital elements in canonical units.  The vector-based
 * constructor and several mutator methods use this function, but
 * programmers should not use it directly in their own code.
 */
void
Traj::randv()
{
    // Set up angles for certain special case orbits.
    // This assumes that if special non-classical orbital
    // Elements are needed to specify the orbit they will
    // have already been defined.

    // Calculate semilatus rectum or semiparameter
    double p = a * ( 1 - e * e );

    // Calculate some temporary values
    double cos_f = cos( f );
    double sin_f = sin( f );
    double temp = p / ( 1.0 + e * cos_f );

    // Calculate position and velocity in PQW frame of reference.
    Vec3 r_pqw( ( temp * cos_f ), ( temp * sin_f ), ( 0.0 ) );

    if ( fabs( p ) < SMALL )
        p = SMALL;

    Vec3 v_pqw( ( -sin_f / sqrt( p ) ), ( ( e + cos_f ) / sqrt( p ) ), 0.0 );

    // Transform orbital frame of reference to geocentric equitorial
    r = rotZ( rotX( rotZ( r_pqw, w ), i ), raan );

    v = rotZ( rotX( rotZ( v_pqw, w ), i ), raan );

    // Fill in e_ h_ and n_vectors.
    double rr = norm( r );

    double vv = norm( v );

    e_vector = ( ( vv * vv - 1.0 / rr ) * r - dot( r, v ) * v ); // CANONICAL UNITS ONLY

    h_vector = cross( r, v );

    n_vector = cross( Vec3( 0, 0, 1 ), h_vector );

    // Everything is solved except for the other anomalies
    // and special-case orbital elements.
    // randv and elorb share this stuff, so it has its own function.
    anomalies();

    special();
} // end randv

/**
 * Calculates classical orbital elements, given radius and velocity
 * in canonical units.  The elements-based contructor and several mutator
 * methods use this function, but programmers should not use it directly
 * in their own code.
 */
void
Traj::elorb( void )
{
    double frac; // temporary variable for bounds checking
    h_vector = cross( r, v );
    n_vector = cross( Vec3( 0, 0, 1 ), h_vector );
    double vv = norm( v );
    double rr = norm( r );
    double hh = norm( h_vector );
    double nn = norm( n_vector );

    // Specific Mechanical Energy, greek letter ksi.  Units: ER^2 / TU^2
    double ksi = vv * vv / 2.0 - 1.0 / rr; // CANONICAL UNITS ONLY
    // double ksi = vv*vv/2.0 - MU/rr // general for non-canonical

    // Semimajor axis.  Units: ER
    a = -1.0 / ( 2.0 * ksi );
    // a = -MU / (2.0 * ksi); // general for non-canonical

    // Eccentricity.  Unitless.
    e_vector = ( ( vv * vv - 1.0 / rr ) * r - dot( r, v ) * v ); // CANONICAL UNITS ONLY
    // Otherwise (S.I. units, MU = 398600 km^3 / s^2)
    // eVector = (1/MU) * ((vv*vv - MU/rr)*r - dot(r,v)*v);
    e = norm( e_vector );

    // Inclination.  Units: Radians
    // Quadrant check is not necessary.
    frac = h_vector.getZ() / hh;

    if (frac > 1.0)
        frac = 1.0;

    if (frac < -1.0)
        frac = -1.0;

    i = acos( frac );

    // RA of ascending node.  Units: Radians.
    frac = n_vector.getX() / nn;

    if (frac > 1.0)
        frac = 1.0;

    if (frac < -1.0)
        frac = -1.0;

    raan = acos( frac );

    // Quadrant check!
    if ( n_vector.getY() < 0 )
        raan = 2 * M_PI - raan;

    // Argument of Perigee. Units: Radians.
    frac = dot( n_vector, e_vector ) / ( nn * e );

    if (frac > 1.0)
        frac = 1.0;

    if (frac < -1.0)
        frac = -1.0;

    w = acos( frac );

    // Quadrant check!
    if ( e_vector.getZ() < 0 )
        w = 2 * M_PI - w;

    // True Anomaly.  Units: Radians.
    frac = dot(e_vector, r) / (e * rr);

    if (frac > 1.0)
        frac = 1.0;

    if (frac < -1.0)
        frac = -1.0;

    f = acos( frac );

    // Quadrant check!
    if ( dot( r, v ) < 0 )
        f = 2 * M_PI - f;

    // Everything is solved except for the other anomalies
    // and special-case orbital elements.
    anomalies();

    special();
}

/**
 * This private function is used by randv() and elorb() to calculate
 * the eccentric and mean anomaly of a trajectory.
 */
void
Traj::anomalies()
{
    if ( e < 1.0 )   // elliptical
    {
        double num = tan( f / 2 );
        double denom = sqrt( ( 1 + e ) / ( 1 - e ) );
        E = 2 * atan2 ( num, denom );

        if ( E < 0.0 )
            E = 2 * M_PI + E;

        M = E - e * sin( E );
    } else if ( fabs( e - 1.0 ) < SMALL )
    { // parabolic
        E = tan( f / 2 ); // hyperbolic anomaly is simple
    }

    else
    { // hyperbolic
        E = 2 * atanh ( sqrt( ( e - 1 ) / ( e + 1 ) ) * tan( f / 2 ) ) ;
        M = e * sinh( E ) - E;
    }

    argLat = w + f;
}

/**
 * This private function is used by randv() and elorb() to handle
 * special case orbits that require non-standard orbital elements.
 * Circular inclined orbits will use argument of latitude.
 * Circular equatorial orbits use true longitude.
 * Noncircular equatorial trajectories use longitude of periapsis.
 */
void
Traj::special()
{
    double rr = norm( r );
    double frac;    // for bounds checking.
    // Argument of Latitude (only for circular inclined orbits)

    if ( e < SMALL )
    {
        frac = dot( n_vector, r ) / ( norm( n_vector ) * rr );

        if (frac > 1.0)
            frac = 1.0;

        if (frac < -1.0)
            frac = -1.0;

        argLat = acos( frac );

        // Quadrant check!
        if ( r.getZ() < 0 )
            argLat = 2 * M_PI - argLat;
    }

    // True Longitude (only for circular equatorial orbits)
    if ( ( e < SMALL ) && ( i < SMALL ) )
    {
        frac = r.getX() / rr;

        if (frac > 1.0)
            frac = 1.0;

        if (frac < -1.0)
            frac = -1.0;

        lonTrue = acos( frac );

        // Quadrant check!
        if ( r.getY() < 0 )
            lonTrue = 2 * M_PI - lonTrue;
    }

    // Longitude of Periapsis (only for non-circular equatorial trajs)
    if ( i < SMALL )
    {
        frac = e_vector.getX() / e;

        if (frac > 1.0)
            frac = 1.0;

        if (frac < -1.0)
            frac = -1.0;

        lonPer = acos( frac );

        // quadrant check
        if ( e_vector.getY() < 0 )
            lonPer = 2 * M_PI - lonPer;
    }
}
