// A simple GNU Scientific Library-based driver to compute the Blasius function
// via implicit integration.  Licensed as GPL3 per GSL dependency.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_minmax.h>
#include <gsl/gsl_odeiv2.h>

// Blasius equation is f''' = - f * f'' / 2 producing a 3D first-order system
enum { ndim = 3 };

// Evaluate the Blasius function RHS
static inline int
forward (double eta, const double f[ndim], double df[ndim], void *params)
{
    (void) eta;
    (void) params;
    df[0] = f[1];
    df[1] = f[2];
    df[2] = - f[0] * f[2] / 2;
    return GSL_SUCCESS;
}

// Use \tilde{\eta} = -\eta to reverse evolution for increasing \tilde{\eta}
static inline int
backward (double eta, const double f[ndim], double df[ndim], void *params)
{
    const int status = forward(-eta, f, df, params);
    for (int i = 0; i < ndim; ++i) df[i] = -df[i];
    return status;
}

// Use GSL ODEIV2 to advance the solution outputting at regular intervals
int
main (int argc, const char *argv[])
{
    // Process optional parameters from command line
    const double etaf   = (argc > 1) ? atof(argv[1]) : 13.6;
    const double deleta = (argc > 2) ? atof(argv[2]) :  0.2;
    const double tol    = (argc > 3) ? atof(argv[3]) : GSL_DBL_EPSILON;

    // Define the problem, tolerance, and initial step size including
    // performing forward or backwards integration from eta=0 initial condition
    const gsl_odeiv2_system sys = {etaf > 0 ? forward : backward, 0L, ndim, 0L};

    // Do /not/ blindly trust the tolerance as Blasius is notoriously difficult
    gsl_odeiv2_driver* const d = gsl_odeiv2_driver_alloc_yp_new (
            &sys, gsl_odeiv2_step_rk8pd, tol*10, tol, 0.0);

    // Initial condition from equation 11 of http://arxiv.org/abs/1006.3888
    // Having a high-precision initial condition avoids shooting approaches
    double eta = 0.0;
    double f[ndim] = { 0.0, 0.0, 0.33205733621519630 };

    // Integrate until either final time achieved or error encountered
    int err = 0;
    printf("%23s  %23s  %23s  %23s\n", "eta", "f", "fp", "fpp");
    const double aetaf = fabs(etaf);
    for (int niter = 0; eta < aetaf; ++niter) { // Odd, but avoids drift
        err = gsl_odeiv2_driver_apply(d, &eta, GSL_MIN(aetaf, niter*deleta), f);
        if (err) {
            fprintf(stderr, "At %g encountered error %d: %s\n",
                    GSL_SIGN(etaf)*eta, err, gsl_strerror(err));
            break;
        }
        printf("%23.16e  %23.16e  %23.16e  %23.16e\n",
               GSL_SIGN(etaf)*eta, f[0], f[1], f[2]);
    }

    gsl_odeiv2_driver_free (d);
    return err;
}
