blasius
=======

The file blasius.c is a simple GNU Scientific Library-based driver to compute
the Blasius function via high order, adaptive integration.  Some color
commentary is available in a post discussing the code:
http://agentzlerich.blogspot.com/2013/10/generating-blasius-boundary-layer.html

For comparison purposes, data taken from a paper by B. D. Ganapol is included
in an accompanying Octave file called blasius.m.  On my system, the output of
the C source run as './blasius 8.8 0.2' matches Ganapol's data at all of his
reported points to at least 8 digits.
