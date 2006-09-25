/*-
* Copyright 2006 (c) Ted Stodgell. All rights reserved.
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
* $Id: Kepler.h,v 1.5 2006/09/25 16:22:20 trs137 Exp $
*
* Contributor(s):  Ted Stodgell <trs137@psu.edu>
*                  David Vallado <valldodl@worldnet.att.net>
*/

#include "Orbgnosis.h"
#include "Stumpff.h"
#include "Traj.h"
#include "Vec3.h"
#include <iostream>
#include <math.h>
using namespace std;

/**
 * Solve Kepler's problem.  Given a state vector (Traj) and a time interval,
 * find the state vector after the time interval has elapsed.  This function
 * calculates r and v vectors and does not consider J2.  If you need J2, 
 * use PKepler() instead.
 * @param traj_0 the initial trajectory at time zero.
 * @param t amount of time, in canonical units.
 *
 * kepler() will throw an integer exception in some cases:
 * If it exceeds the iteration limit it throws 1.
 * If it converges, but the F&G transformation is out of tolerance, it throws
 * a 2.
 * Obviously, programmers are encouraged to wrap calls to kepler()
 * in a try-catch block.
 */
Traj
kepler ( Traj traj_0, double t )
{
    //srand (time (NULL));
    if (t < 0) 
    {
        cout << "Kepler needs time > 0." << endl;
        exit(1);
    }
    if ( fabs( t ) <= SMALL )
    {
        cout << "Kepler: time was zero.  No movement." << endl;
        return traj_0; // Zero time, so no movement.
    } else {
        // set up local variables
        double r0, v0;  // initial radius and velocity (magnitudes only)
        Vec3 rfinal, vfinal;      // final radius and velocity (vectors)
        double F, G;    // universal variable f and g expressions
        double Fdot, Gdot; // F and G's rates of change
        double Xold;    // universal variable
        double Xnew;    // universal variable
        double Xold2;   // Xold squared
        double Xnew2;   // Xnew squared
        double Znew;    // new value of Z
        double C2new;   // Stumpff C2 value
        double C3new;   // Stumpff C3 value
        double tnew;    // new time
        double rdotv;   // r0 dot v0;
        double a;       // semimajor axis;
        double alpha;   // 1 / (semimajor axis)
        double ksi;     // specific mechanical energy
        double period;  // orbital period
        double S, W;    // variables for parabolic special case
        double Rval;
        double temp;
        double adjust;  // variable step adjuster.
        int counter = 0;
        int adj_ctr = 0;
        const int limit = 400;  // iteration limit (default = 40)
        Traj result;

        r0 = norm( traj_0.get_r() ); // current radius
        v0 = norm( traj_0.get_v() ); // current velocity
        Xold = 0.0;
        Znew = 0.0;
        rdotv = dot( traj_0.get_r(), traj_0.get_v() );

        ksi = ( 0.5 * v0 * v0 ) - ( 1.0 / r0 );     // canonical!
        alpha = -2.0 * ksi;
        a = traj_0.get_a();

        // Set up initial guess for Xold.
        if ( alpha >= 0.0001 ) // was (alpha >= SMALL)
        {
            //cout << "**** Kepler is working on an ellipse ****" << endl;
            period = 2 * M_PI * sqrt( pow( fabs( a ), 3.0 ) );
            if ( fabs( t ) > fabs( period ) )
            {
                // cout << "Kepler is dealing with >1 revolution." << endl;
                t = fmod ( t, period ); // multirev
            }
            if ( fabs( alpha - 1.0 ) > 0.5 )
                Xold = t * alpha;
            else
                Xold = t * alpha * 0.97; // first guess can't be too close
        } else {
            if ( fabs(alpha) < 0.0001 ) // was (fabs(alpha) < SMALL )
            {
                // Parabola
                //cout << "**** Kepler is working on a parabola ****" << endl;
                double h = norm( traj_0.get_h_vector() );
                double p = h * h;
                S = 0.5 * ( M_PI / 2.0 - atan( 3.0 * sqrt( 1.0 / ( p * p * p ) ) * t ));
                W = atan( pow( tan(S), 1.0 / 3.0 ) );
                Xold = sqrt(p) * ( 2.0 * ( 1.0 / tan( 2.0 * W ) ) );
                // alpha = 0.0;  // experiment!!!
            } else {
                // Hyperbola
                // This only works correctly for positive t.
                //cout << "**** Kepler is working on a hyperbola ****" << endl;
                temp = -2.0 * t /
                       ( a * ( rdotv + sqrt( -a ) * ( 1.0 - r0 * alpha ) ) );
                Xold = sqrt( -a ) * log( temp );
            }
        }

        while ( 1 )    // XXX ugly loop
        {
            Xold2 = Xold * Xold;
            Znew = Xold2 * alpha;
            C2new = stumpff_C2( Znew );
            C3new = stumpff_C3( Znew );

            tnew = Xold2 * Xold * C3new + rdotv * Xold2 * C2new +
                   r0 * Xold * ( 1.0 - Znew * C3new );

            Rval = Xold2 * C2new + rdotv * Xold * ( 1.0 - Znew * C3new ) +
                   r0 * ( 1.0 - Znew * C2new );

            if ((a > 0.0) &&
                (fabs(Xnew) > 2*M_PI*sqrt(a)) &&
                (ksi < 0.0))
            {
                /* XXX
                 * This adjusts the step size if things are going badly.
                 * Rval is multiplied by a random number between 7 and 10.
                 * Why random?  Because I can't think of a good constant value.
                 * Vallado recommends 7 to 10.
                 */
                adjust = 3.0*((double)rand()/((double)(RAND_MAX)+(double)(1)))+7.0;
                Xnew = Xold + (t - tnew) / (Rval * adjust);
                adj_ctr++;
            } else {
                Xnew = Xold + ( t - tnew ) / Rval;
            }
            counter++;
            //cout << tnew-t << endl;
            Xold = Xnew;
            if (( fabs( tnew - t ) < SMALL ) || ( counter >= limit )) break;
        }  // end while

        // Update Znew, C2new and C3new!
        // Vallado's original code doesn't have this.
        // Not doing this will cause small errors, 
        // especially with parabolic cases.
        Xold2 = Xold * Xold;
        Znew = Xold2 * alpha;
        C2new = stumpff_C2( Znew );
        C3new = stumpff_C3( Znew );
        //---------------------------------------------------

        //cout << "Kepler finished iterating after " << counter << " times." << endl;
        // Calculate position and velocity vectors at new time
        Xnew2 = Xnew * Xnew;
        F = 1.0 - ( Xnew2 * C2new / r0 );
        G = t - Xnew2 * Xnew * C3new;
        rfinal = F * traj_0.get_r() + G * traj_0.get_v();
        Gdot = 1.0 - ( Xnew2 * C2new / norm( rfinal ) );
        Fdot = ( Xnew / ( r0 * norm( rfinal ) ) ) * ( Znew * C3new - 1.0 );
        vfinal = Fdot * traj_0.get_r() + Gdot * traj_0.get_v();
        temp = F * Gdot - Fdot * G;

        if ( counter >= limit ) throw(1);
        if ( fabs( temp - 1.0 ) > 0.00001 ) throw(2);

        // Traj constructor automatically takes care of classical elements.
        result = Traj (rfinal, vfinal);
        result.do_J2_regression(t);
        return result;
    }
}
