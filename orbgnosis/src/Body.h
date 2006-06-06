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
 * $Id: Body.h,v 1.7 2006/06/06 15:26:16 trs137 Exp $
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 */

#ifndef _BODY_H_
#define _BODY_H_
#include "Vector.h"
#include "Traj.h"

/**
 * The Body class represents any massive body in space.
 * Mass properties and state vector are stored here.
 */
class Body
{
    private:
        double  mass;
        Vector  position;
        Vector  velocity;
        Vector  ang_vel;        // angular velocity, rad/s
        Vector  moments;        // principle moments of inertia
        Traj    trajectory;     // current trajectory of the Body.

    public:
                Body (void);    // Default ctor
                Body (const Body&);       // copy constructor
                Body& operator = (Body);  // copy assignment operator
        virtual ~Body (void);             // destructor

        double  kineticEnergy (void);     // kinetic energy of body
        double  transMomentum (void);     // translational momentum
        Vector  angMomentum (void);       // angular momentum

        void    move (Vector);
        void    accelerate (Vector);
        Vector  findForce (void);         // sum of all forces on Body
};

#endif /* _BODY_H_ */
