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
 * $Id: Body.cpp,v 1.6 2006/06/06 15:07:15 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 */

#include "Body.h"

/**
 * Default Body constructor with no args.
 * Sets name to "NO NAME", mass to 1.0, and everything else to zero.
 */
Body::Body (void)
{
    name = "NO NAME                       ";
    mass = 1.0;
    position = Vector(0.0, 0.0, 0.0);
    velocity = Vector(0.0, 0.0, 0.0);
    ang_vel  = Vector(0.0, 0.0, 0.0);
    moments  = Vector(0.0, 0.0, 0.0);
    //trajectory = Traj(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}
/**
 * The Body destructor.
 */
Body::~Body (void)
{
    //
}

/**
 * The Body copy constructor.
 */
Body::Body (const Body& copy)
{
    mass = copy.mass;
    position = copy.position;
    velocity = copy.velocity;
    ang_vel = copy.ang_vel;
    moments = copy.moments;
    //trajectory = copy.trajectory;
}

/**
 * The Body copy assignment operator.
 */
Body&
Body::operator = (Body b)
{
    name = b.name;
    mass = b.mass;
    position = b.position;
    velocity = b.velocity;
    moments = b.moments;
    trajectory = b.trajectory;
    return *this;
}

/**
 * Returns the kinetic energy of the Body.
 */
double
Body::kineticEnergy (void)
{
    double speed = norm(velocity);
    return (0.5 * mass * speed * speed );
}

/**
 * Returns the translational momentum of the Body.
 */
double
Body::transMomentum (void)
{
    double speed = norm(velocity);
    return (mass * speed);
}

/**
 * Returns the angular momentum Vector of the Body.
 */
Vector
Body::angMomentum (void)
{
    return (cross(ang_vel, moments));
}

/**
 * Moves the Body by vector addition.
 * @param vin is the Vector to be added to the Body's position.
 */
void
Body::move (Vector vin)
{
    position = position + vin;
}

/**
 * Accelerates the Body by vector addition.
 * @param vin is the Vector to be added to the Body's velocity.
 */
void
Body::accelerate (Vector vin)
{
    velocity = velocity + vin;
}

/**
 * Returns the sum of all forces acting upon the Body.
 */
Vector
Body::findForce (void)
{
    return Vector(0,0,0); //TODO
}
