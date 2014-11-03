This directory contains source files for
the variations on the deflation scheme.

Files
======
KEY: (T) - essentially completely broken code
(B) - mostly working, but buggy and breaks sometimes
(I) - works, but doesn't do what you think it does
(W) - works

hessen.m (W):
  Uses householder transformations to
  transform a matrix into hessenberg form.
  Can do it both upwards and downwards and catches hermitian matrices.

pltMat.m (W):
  Generates a plot of a matrix

trainBust.m (W):
  Chases two bulge trains towards center from opposite ends, then explodes them.
  Will also chase bulges to the top and bottom

schurAd.m (W):
  Computes the derivative of both S and U from a Schur decomposition
  (A.U = U.S)

bottom/adEvals.m (T):
  Calculates eigenvalues using aggresive early deflation.
  NOTE - this doesn't perform as advertised, currently has
  some logic built in but far from working/finished

bottom/agErlyDef.m (W):
  Chases a train of bulges to the bottom, creates a spike,
  then cleans up the spike.
  The shifts used are not the traditional qr shifts.

testScripts/ (W?):
  functions to perform tests on a variety of different functions and situations

Notes on Behavior
=================

trainBust:
  If the shifts used are from eig(A), the following occurs in a single pass:
  -pushing 'against the grain' - yeilds ~ as many zeros in spikes as there are shifts
  -pushing 'with the grain' - yeilds a variable amount of zeros in spikes, near 0 with extremal eigenvalues
  -pushing bulges towards middle - yields very few zeros in spikes
  'with the grain' -> both the hessenberg formation step and the bulges were pushed in the same direction
  'against the grain' -> the hessenberg formation step and the bulges were pushed opposite directions
