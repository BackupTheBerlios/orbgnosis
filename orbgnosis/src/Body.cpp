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
 * $Id: Body.cpp,v 1.4 2006/03/19 22:05:34 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 */

#include "Body.h"

Body::Body (void)
{
    name = "NO NAME                       ";
    mass = 1.0
    position = new Vector(); // 0, 0, 0
    velocity = new Vector(); // 0, 0, 0
}

Body::~Body (void)
{
    //
}

Body::Body (const Body& copy)
{
    name = copy.name;  // this is wrong
    mass = copy.mass;
    position = new Vector();
    position = copy.position;
    velocity = new Vector();
    velocity = copy.velocity;
}

double
Body::kineticEnergy (void)
{
    double speed = velocity.norm();
    return (0.5 * mass * speed * speed );
}

double
Body::transMomentum (void)
{
    double speed = velocity.norm();
    return (mass * speed);
}

void
Body::move (Vector vin)
{
    position.add(vin);
}

void
Body::accelerate (Vector vin)
{
    velocity.add(vin);
}

Vector
Body::findForce (void)
{
    //TODO
}
