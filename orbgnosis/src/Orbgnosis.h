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
* $Id: Orbgnosis.h,v 1.21 2006/09/12 18:25:16 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

/** @file
 * This header file defines constants and has some utitlity functions.
 */

#ifndef _ORBGNOSIS_H_
#define _ORBGNOSIS_H_

#ifdef INF
#undef INF
#endif
#define INF 1.0e14    //!< A very big number

#define EPS 1.0e-14   //!< A very small number, used in NSGA2 
//#define E  2.71828182845905  //!< The natural number e
//#define PI 3.14159265358979   use M_PI from <math.h> instead
#define GNUPLOT_COMMAND "gnuplot -persist" //!<  Used for drawing with gnuplot

#define SMALL 1.0e-8   //!< Tolerance and general-purpose "really small number"
#define MU 398600.4418 //!< the gravitational parameter of Earth in km<sup>3</sup> s<sup>-2</sup> in S.I. units.
#define ROOTMU 631.3481     //!< the square root of MU in S.I. units.

/* A NOTE ON "CANONICAL" TIME
 * Time and length units are chosen so that MU is unity.  This
 * minimizes roundoff and truncation error, as well as saves floating
 * point operations involving MU.
 *
 * In geocentric trajectories, Earth's radius is a convenient unit of
 * length, so that's what we will use.
 *
 * Now, suppose there is a circular orbit with radius = 1 length unit.
 * Its semimajor axis is also 1 length unit.
 * Its orbital period is 2 * pi * sqrt( a*a*a / MU ), and since we wish MU
 * to be one, we may solve the expression:
 * In canonical time units, a trajectory with a semimajor axis of 1 length 
 * unit has a period of exactly two pi time units.  Knowing this, we may
 * easily find the length of a canonical time unit in terms of seconds
 * or minutes.
 */
#define ER 6378.137    //!< Earth radius in km... (km per canonical length unit)
#define TU_SEC 806.81112382429    //!< seconds per canonical time unit
#define TU_MIN 13.44685206374     //!< minutes per canonical time unit

#define nan std::numeric_limits<double>::quiet_nan(); //!< ifndef NAN...

#endif /* _ORBGNOSIS_H_ */
