#does automatic differentiation with newton's method to minimize
#subdiagonals. Jacobian tends to become singular pretty quickly with random
#matrices
function [itTilSplit] = midDefAD(H,s)
  MAXITER = 150;
  itTilSplit = 0;
  n = length(H);
  n1 = n;
  tol = sqrt(eps);
  shifts = eig(H(n-s+1:n,n-s+1:n));
  shifts = 10*randn(s,1);
  seed = eye(s);
  fprintf('starting msqriAD loop\n');
  mids = ceil(s/2)-2:ceil(s/2)+s+2;
  do
    
    n = length(H);
    itTilSplit++;
    fprintf('%d ',itTilSplit);
    [B, ~,Jac,f,INFO] = trainBust(H,real(shifts),seed,false,false,'m',true);
    n = length(H);
    f = f';
    Jac = Jac';
    #only consider the s middle shifts
    f = f(mids);
    Jac = Jac(mids,:);
    if(iscell(B))
      return;
    end
  
    #Get low rank stuff
    [U,S,V] = svd(Jac,1);
    S = diag(S);
    S(S < .1*S(1,1)) = [];
    U = U(:,1:length(S));
    V = V(:,1:length(S));

    up = V*diag(1./S)*U'*f;
    if(norm(up)/norm(shifts) < 1e-4)
      shifts = randn(length(shifts),1);
      pltMat(B);
    else
      #Newton
      shifts = shifts - V*diag(1./S)*U'*f
    end

  until(itTilSplit > MAXITER)
end
