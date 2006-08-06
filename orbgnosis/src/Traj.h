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
* $Id: Traj.h,v 1.9 2006/08/06 00:33:49 trs137 Exp $
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
        Traj (void); // defaults to all zeroes

        // This version of the constructor will call randv();
        Traj (double,    // a
              double,    // e
              double,    // i
              double,    // raan
              double,    // w
              double);  // f

        // This version of the constructor will call elorb().
        Traj (Vec3,   // r
              Vec3);  // v

        virtual ~Traj (void);

        Traj (const Traj&); // copy constructor

        Traj& operator = (Traj); // copy assignment operator

        void print (void);

        double get_a (void);
        double get_e (void);
        double get_i (void);
        double get_raan (void);
        double get_w (void);
        double get_f (void);
        Vec3 get_r (void);
        Vec3 get_v (void);

        void set_a (double);
        void set_e (double);
        void set_i (double);
        void set_raan (double);
        void set_w (double);
        void set_f (double);
        void set_r (Vec3);
        void set_v (Vec3);

        void randv (void); // Calculates r and v vectors from classical elements.
        void elorb (void); // Calculates classical elements from 2 vectors.

    private:
        double a;       //!< Semimajor Axis (length).
        double e;       //!< Eccentricity (dimensionless).
        double i;       //!< Inclination  (radians).
        double raan;    //!< Right Ascension of the Ascending Node (radians).
        double w;       //!< Argument of Perigee (radians).
        double f;       //!< True Anomaly (radians).
        Vec3 r;         //!< Radius (ER).
        Vec3 v;         //!< Velocity (ER/TU).
};

#endif /* _TRAJ_H_ */
