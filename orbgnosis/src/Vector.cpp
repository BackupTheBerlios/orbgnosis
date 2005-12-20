/*
 * Copyright 2005 Ted Stodgell. All rights reserved.
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
 * $Id:
 *
 * Contributor(s):  Ted Stodgell <trs137@psu.edu>
 */

#include "Vector.h"
#include <math.h>

Vector::Vector (void)
{
    x = y = z = 0;
}

Vector::Vector (double xin, double yin, double zin)
{
    x = xin;
    y = yin;
    z = zin;
}

Vector::~Vector (void)
{
    //
}

void
Vector::setX (double xin)
{
    x = xin;
}

void
Vector::setY (double yin)
{
    y = yin;
}
void
Vector::setZ (double zin)
{
    z = zin;
}

void
Vector::add (Vector vin)
{
    x = x + vin.x;
    y = y + vin.y;
    z = z + vin.z;
}

double
Vector::getX (void)
{
    return x;
}

double
Vector::getY (void)
{
    return y;
}

double
Vector::getZ (void)
{
    return z;
}

double
Vector::norm (void)
{
    return ( sqrt (x*x + y*y + z*z) );
}
