#Compute the derivative of the schur decomposition
function [Qd] = schurAd(Ad, Q, T,tol = 1e-12)
  B = Q' * Ad * Q;
  P = recurSolve(zeros(size(T)),0,T,B);
  P = P - P';
  #compute the Q derivative
  Qd = Q*P;
  Qd(abs(Qd) < tol) = 0;#eliminate floating point zeros
  
  #compute the T derivative
  #Td = T*P - P*T + B
endfunction

#calculates the zero on the subdiagonal
#which is closest to the center
function [n] = spltPoint(A,tol=1e-12)
  idx = find(abs(diag(A,-1)) <= tol);#find all zeros
  if(isempty(idx))#if no split point
    n = -length(A);
    return
  endif
  [~,pos] = sort(abs(idx - length(A)/2));#sort them by closest to middle, and return index of them
  n = idx(pos(1));
  return;
endfunction

#Solves the system recursively
function [P] = recurSolve(P,n,T,B)
  m = spltPoint(T);
  #if finished with computation
  if( m == 0)
    return;
  elseif( m == -2)
    P(n + 2, n + 1) = .5 * (B(2,2) - B(1,1)) / (T(1,2) + T(2,1));
    return;
  endif

  #splitting indices
  t = 1:m;
  b = m+1:length(T);

  X = syl(T(b,b), -T(t,t), B(b,t));
  [n1 n2] = size(X);
  P( n+m + (1:n1) , n + (1:n2) ) = X;

  #recurse on remaining blocks
  P = recurSolve(P, n, T(t,t), B(t,t) + B(t,b)*X);
  P = recurSolve(P,n+m,T(b,b), B(b,b) - X*T(t,b));
  return;
endfunction

