errBlock = 0;
errDiag = 0;

for i=1:30
  #This tests whether the derivative given is actually the derivative of the schur decomposition.
  A = randn(50);
  #seed the derivative of A as random
  for j=1:10
    Ad = randn(50);

    #form the schur decomposition
    [Q,T] = schur(A);
    [Qd,Td] = schurAd(Ad,Q,T);

    #figure out how far off it is
    errBlock += norm(Ad*Q + A*Qd - Qd*T - Q*Td,'fro');
 
    #form the schur decomposition
    [Q,T] = schur(A,'complex');#gaurantee actually diagonal
    [Qd,Td] = schurAd(Ad,Q,T);

    #figure out how far off it is
    errDiag += norm(Ad*Q + A*Qd - Qd*T - Q*Td,'fro');
  end
  i
end
errBlock /= 300
errDiag /= 300
