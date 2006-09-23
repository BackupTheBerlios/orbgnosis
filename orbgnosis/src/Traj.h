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
* $Id: Traj.h,v 1.24 2006/09/23 04:03:44 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

#ifndef _TRAJ_H_
#define _TRAJ_H_
#include "Vec3.h"

/**
 * The extrinsic properties of a body: its orbit or trajectory through space.
 * The Traj class represents a geocentric trajectory and
 * contains classical Keplerian orbital elements needed to specify that
 * trajectory uniquely.  Target satellites are defined in part by their
 * trajectories.  Transfer arcs between targets are also trajectories.
 */

class Traj
{

    public:
        Traj ( void ); // defaults to all zeroes

        // This ctor calls randv() to fill in the missing data.
        Traj ( double,      // a
               double,      // e
               double,      // i
               double,      // raan
               double,      // w
               double );   // f

        // This ctor calls elorb() to fill in the missing elements.
        Traj ( Vec3,     // r
               Vec3 );  // v

        virtual ~Traj ( void );       // dtor
        Traj ( const Traj& );         // copy constructor
        Traj& operator = ( Traj );    // copy assignment

        void print ( void );  // Prints trajectory parameters to stdout.
        void print_El (void); // Prints just the 6 classical elements one line.

        // Accessors
        double get_a ( void );
        double get_e ( void );
        double get_i ( void );
        double get_raan ( void );
        double get_w ( void );
        double get_f ( void );
        Vec3 get_r ( void );
        Vec3 get_v ( void );
        double get_E ( void );
        double get_M ( void );
        double get_argLat ( void );
        double get_lonTrue ( void );
        double get_lonPer ( void );
        Vec3 get_e_vector ( void );
        Vec3 get_h_vector ( void );
        Vec3 get_n_vector ( void );

        // Mutators
        void set_elorb ( double, double, double, double, double, double );
        void set_randv ( Vec3, Vec3 );
        void set_a ( double );
        void set_e ( double );
        void set_i ( double );
        void set_raan ( double );
        void set_w ( double );
        void set_f ( double );
        void set_r ( Vec3 );
        void set_v ( Vec3 );
        void set_M ( double );

    private:
        // SIX CLASSICAL ORBITAL ELEMENTS:
        double a;       //!< Semimajor Axis (length).
        double e;       //!< Eccentricity (dimensionless).
        double i;       //!< Inclination  (radians).
        double raan;    //!< Right Ascension of the Ascending Node (radians).
        double w;       //!< Argument of Perigee (radians).
        double f;       //!< True Anomaly (radians).

        // GEOCENTRIC EQUATORIAL STATE VECTOR:
        Vec3 r;         //!< Radius (ER).
        Vec3 v;         //!< Velocity (ER/TU).

        // OTHER ORBITAL ELEMENTS:
        double E;        //!< Eccentric, Parabolic, or Hyperbolic anomaly (radians).
        double M;        //!< Mean Anomaly (radians).
        double argLat;   //!< Argument of Latitude (radians).
        double lonTrue;  //!< True Longitude (radians).
        double lonPer;   //!< Longitude of Periapsis (radians).

        // OTHER VECTORS:
        Vec3 e_vector;  //!< Eccentricity
        Vec3 h_vector;  //!< Specific angular momentum
        Vec3 n_vector;  //!< Node vector

        // Private methods.
        void randv ( void );     // Calculates r and v vectors from classical elements.
        void elorb ( void );     // Calculates classical elements from 2 vectors.
        void anomalies ( void ); // Routines common to randv() and elorb().
        void special ( void );   // Calculates orb elements for special case orbits.
};

#endif /* _TRAJ_H_ */
