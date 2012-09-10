GraviDy
=======

A quasi-keplerian n-body gravitational system integrator.
----------------------------------------------------------

The astrophysics "N-body problem" consist in a dynamical
N celestial bodies model which interact gravitationally.

Solving this problem implies understand and predict
the movement of the system components, for example,
the Solar System, a stars cumulus, or a galaxy.

To systems with more than 3 bodies,
which are most of the interesting astronomical systems,
the problem do not have an analytic solution,
and it is necessary to use computational methods,
which are progresively more expensives to largest values of N.


This projects pretends to use the GPU technology to solve
the N-body problem of a "quasi-keplerian" system,
in other words, a system in which their dynamics is dominated
by a massive central object, but the interactions between
the others bodies are no totally negligible.
