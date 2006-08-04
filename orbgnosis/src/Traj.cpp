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
* $Id: Traj.cpp,v 1.8 2006/08/04 02:49:11 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

#include <iostream>
#include "Orbgnosis.h"
#include "Traj.h"
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
        f(0.0)
{
    //std::cout << "Traj constructor called with zero args.\n";
}

/**
 * The Traj constructor.
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
        f(fin)
{
    //std::cout << "Traj constructor called with 6 args.\n";
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
        f(copy.f)
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
void Traj::set_a (double ain) { a = ain; }
void Traj::set_e (double ein) { a = ein; }
void Traj::set_i (double iin) { a = iin; }
void Traj::set_raan (double raanin) { a = raanin; }
void Traj::set_w (double win) { a = win; }
void Traj::set_f (double fin) { a = fin; }
