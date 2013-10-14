// A simple GNU Scientific Library-based driver to compute the Blasius function
// via implicit integration.  Licensed as GPL3 per GSL dependency.

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_odeiv2.h>

// Blasius equation is f''' = - f * f'' / 2 producing a 3D first-order system
enum { ndim = 3 };

// Evaluate the Blasius RHS
static inline int
func (double eta, const double f[], double df[], void *params)
{
    (void) eta;
    (void) params;
    df[0] = f[1];
    df[1] = f[2];
    df[2] = - f[0] * f[2] / 2;
    return GSL_SUCCESS;
}

// Evaluate the Blasius Jacobian and delegate for the RHS
// http://www.gnu.org/software/gsl/manual/html_node/Defining-the-ODE-System.html
static int
jac (double eta, const double f[], double* d2f, double df[], void *params)
{
    d2f[/*0*/  + 0] = 0;
    d2f[/*0*/  + 1] = 1;
    d2f[/*0*/  + 2] = 0;
    d2f[  ndim + 0] = 0;
    d2f[  ndim + 1] = 0;
    d2f[  ndim + 2] = 1;
    d2f[2*ndim + 0] = -f[0] / 2;
    d2f[2*ndim + 1] = 0;
    d2f[2*ndim + 2] = -f[2] / 2;

    return func(eta, f, df, params);
}

// Use GSL ODEIV2 to implicitly advance the solution outputting regularly
int
main (void)
{
    // Define the problem, tolerance, and initial step size
    // Do not trust in the tolerance as Blasius is notoriously difficult
    const double tol = GSL_DBL_EPSILON;
    const gsl_odeiv2_system sys = {func, jac, ndim, NULL};
    gsl_odeiv2_driver* const d = gsl_odeiv2_driver_alloc_yp_new (
            &sys, gsl_odeiv2_step_rk4imp, sqrt(tol), tol, 0.0);

    // Prepare output frequency, initial eta, and final eta
    const double deleta = 0.2;
    double       eta    = 0.0;
    double       etaf   = 8.8;

    // Initial condition from equation 11 of http://arxiv.org/abs/1006.3888
    double f[ndim] = { 0.0, 0.0, 0.33205733621519630 };

    int status;
    printf("%22s  %22s  %22s  %22s\n", "eta", "f", "fp", "fpp");
    printf("%22.16e  %22.16e  %22.16e  %22.16e\n", eta, f[0], f[1], f[2]);
    while (eta < etaf)
    {
        status = gsl_odeiv2_driver_apply (d, &eta, eta + deleta, f);
        if (GSL_SUCCESS != status) {
            fprintf(stderr, "At %g encountered error %d: %s\n",
                    eta, status, gsl_strerror(status));
            break;
        }
        printf("%22.16e  %22.16e  %22.16e  %22.16e\n", eta, f[0], f[1], f[2]);
    }

    gsl_odeiv2_driver_free (d);
    return 0;
}
