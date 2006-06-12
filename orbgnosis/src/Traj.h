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
* $Id: Traj.h,v 1.7 2006/06/12 21:22:17 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

#ifndef _TRAJ_H_
#define _TRAJ_H_

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
        Traj (void);
        Traj (double,    // a
              double,    // e
              double,    // i
              double,    // raan
              double,    // w
              double);  // f

        virtual ~Traj (void);

        Traj (const Traj&); // copy constructor

        Traj& operator = (Traj); // copy assignment operator

        double getfoo (void);
        void setFoo (double);

    protected:
        double a;       //!< Semimajor Axis (length).
        double e;       //!< Eccentricity (dimensionless).
        double i;       //!< Inclination  (radians).
        double raan;    //!< Right Ascension of the Ascending Node (radians).
        double w;       //!< Argument of Perigee (radians).
        double f;       //!< True Anomaly (radians).
};

#endif /* _TRAJ_H_ */
