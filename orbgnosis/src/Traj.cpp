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
 * $Id: Traj.cpp,v 1.1 2006/06/05 18:38:58 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 */

#include "Orbgnosis.h"
#include <math.h>

/*
 * This pragma disables 
 * "remark #981: operands are evaluated in unspecified order"
 * on the Intel C/C++ compiler.  ICC warns about this in unneccessary
 * cases.  Leave NO_ICC_981 undefined to get the warnings.
 */
#ifdef NO_ICC_981
#pragma warning (disable:981)
#endif

/**
 * The Traj constructor.
 * All six classical Keplerian elements must be
 * specified in the constructor.
 * @param ain the semimajor axis (length units).
 * @param ein the eccentricity (dimensionless).
 * @param iin the inclination (radians).
 * @param raanin the right ascention of the ascending node (radians)
 * @param win the argument of periapsis (radians).
 * @param fin the true anomaly (radians).
 */
Traj::Traj (double )
{
    a = ain;
    e = ein;
    i = iin;
    raan = raanin;
    w = win;
    f = fin;
}

/**
 * The Traj destructor.
 */
Traj::~Traj (void)
{
    // cout << "Traj destructor called.\n";
}

/**
 * The Traj copy constructor.
 * @param the reference to the Traj to be copied.
 */
Traj::Traj (const Traj& copy)
{
    a = copy.a;
    e = copy.e;
    i = copy.i;
    raan = copy.raan;
    w = copy.win;
    f = copy.f;
}

/**
 * The Tracj copy assignment operator.
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


// Re-enable ICC remark #981
#ifdef NO_ICC_981
#pragma warning (default:981)
#endif
