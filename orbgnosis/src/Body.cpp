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
* $Id: Body.cpp,v 1.11 2006/06/12 21:22:17 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*/

#include "Orbgnosis.h"
#include "Body.h"

/**
 * Default Body constructor with no args.
 * Sets mass and moments to 1.0, and everything else to zero.
 */
Body::Body (void)
    : mass(1.0),
      position(0.0, 0.0, 0.0),
      velocity(0.0, 0.0, 0.0),
      ang_vel(0.0, 0.0, 0.0),
      moments(1.0, 1.0, 1.0)
{
    //cout << "Body constructor called with no args.\n";
}

/**
 * The Body destructor.
 */
Body::~Body (void)
{
    // cout << "Body destructor called\n";
}

/**
 * The Body copy constructor.
 */
Body::Body (const Body& copy)
    : mass(copy.mass),
      position(copy.position),
      velocity(copy.velocity),
      ang_vel(copy.ang_vel),
      moments(copy.moments)
{
    //cout << "Body copy constructor called\n";
}

/**
 * The Body copy assignment operator.
 */
Body&
Body::operator = (Body b)
{
    mass = b.mass;
    position = b.position;
    velocity = b.velocity;
    moments = b.moments;
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
    return Vector(0, 0, 0); //TODO not really necessary yet
}
