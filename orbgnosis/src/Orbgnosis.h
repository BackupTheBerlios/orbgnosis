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
* $Id: Orbgnosis.h,v 1.14 2006/08/07 02:32:26 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

/** @file
 * This header file defines global constants and macros.
 * Any constants not defined in standard headers should be defined here.
 */

#ifndef _ORBGNOSIS_H_
#define _ORBGNOSIS_H_

#ifdef INF
#undef INF
#endif // ifdef INF
#define INF 1.0e14    //!< A big number.

#define EPS 1.0e-14   //!< A very small number, used in NSGA2
//#define E  2.71828182845905  //!< The natural number e
//#define PI 3.14159265358979   use M_PI from <math.h> instead
#define GNUPLOT_COMMAND "gnuplot -persist" //!<  Command for drawing with gnuplot

#define SMALL 1.0e-8   //!< Tolerance and general-purpose "small number".
#define VERY_SMALL 1.0e-16 //!< Tolerance and general-purpse "very small number".
#define ER 6378.137    //!< One canonical length unit.  Earth radius in km.
#define MU 398600.4418 //!< The NON-CANONICAL gravitational parameter of Earth in km<sup>3</sup> s<sup>-2</sup>.  Note that MU = 1 in canonical units.
#define ROOTMU 631.3481     //!< the square root of Mu in SI units..
#define TU_SEC 806.81112382429    //!< One canonical time unit in seconds.
#define TU_MIN 13.44685206374     //!< One canonical time unit in minutes.

#endif /* _ORBGNOSIS_H_ */
