blasius
=======

The file blasius.c is a simple GNU Scientific Library-based driver to compute
the Blasius function via implicit, adaptive integration.

For comparison purposes, data taken from a paper by B. D. Ganapol is included
in an accompanying Octave file called blasius.m.  On my system, the output of
the C source run as './blasius 8.8 0.2' matches Ganapol's data at eta = 8.8 for
all but the final digit he reports in the second derivative of f.
