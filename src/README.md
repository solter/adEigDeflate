This directory contains source files for
the variations on the deflation scheme.

Files
======
Note: most of the header information in the .m files
is outdated and incorrect. This is due to many of
these files essentially being branches of other
files. The comments throughout the code should
be much better maintained and the variable names
are for the most part descriptive.
Cleanup is necesary and recombination
would be a good thing. Furthermore, in the future
please do this with git branching, not creating
many files (my bad for the current state, sorry).

KEY: (T) - essentially completely broken code
(B) - mostly working, but buggy and breaks sometimes
(I) - works, but doesn't do what you think it does
(W) - works

hessen.m (W):
  Uses householder transformations to
  transform a matrix into hessenberg form

pltMat.m (W):
  Generates a plot of a matrix

bottom/adEvals.m (T):
  Calculates eigenvalues using aggresive early deflation.
  NOTE - this doesn't perform as advertised, currently has
  some logic built in but far from working/finished

bottom/agErlyDef.m (I):
  Chases a train of bulges to the bottom, creates a spike,
  then cleans up the spike.

middle/colAndBust.m (I):
  Chases two bulges towards center from opposite ends, then explodes them.
  This doesn't push trains. The shifts aren't the traditional shifts

middle/trainBust.m (W?):
  Chases two bulge trains towards center from opposite ends, then explodes them.
  Will also chase bulges to the top and bottom

middle/eigBySplit.m (T):
  Uses colAndBust to blow up matrix, then extracts the spikes.
  Is supposed to then use this info to update shift strategy and
  keep going until it can split the matrix, but currently is
  missing this logic.

Intended function usage
=======================
KEY:
*name - basic function (level 0 function)
->name - level 1 - function utilizing basic function to perform more advanced calculation
->->name - level 2 -  function utilizing level 1 function to perform more advanced calculation

* pltMat :
  for visualization

* hessen :
  for reference. Octave's built in hess function is faster.

* agErlyDef :
  basic train chasing with option of spike

-> adEvals :
  run through agErlyDef a bunch to perform complete eval calculation

* colAndBust:
  see above

-> eigBySplit :
  run through colAndBust to complete eval calculation

