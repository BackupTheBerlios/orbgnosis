/*-
 * Copyright (c) 2005 Ted Stodgell. All rights reserved.
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
 * $Id: Sat.cpp,v 1.1 2006/06/06 19:56:35 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 */

#include "Orbgnosis.h"
#include "Sat.h"

/**
 * Default Sat constructor with no args.
 * Sets mass to 1.0, and everything else to zero.
 */
Sat::Sat (void)
{
    // Nothing in Sat yet.
}

/**
 * The Sat destructor.
 */
Sat::~Sat (void)
{
    // cout << "Sat destructor called\n";
}

/**
 * The Sat copy constructor.
 */
Sat::Sat (const Sat& copy)
{
    // Body
    mass = copy.mass;
    position = copy.position;
    velocity = copy.velocity;
    ang_vel = copy.ang_vel;
    moments = copy.moments;
    // Traj
    a = copy.a;
    e = copy.e;
    i = copy.i;
    raan = copy.raan;
    w = copy.w;
    f = copy.f;
}

/**
 * The Sat copy assignment operator.
 */
Sat&
Sat::operator = (Sat b)
{
    mass = b.mass;
    position = b.position;
    velocity = b.velocity;
    moments = b.moments;
    return *this;
}
