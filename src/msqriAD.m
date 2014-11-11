#does automatic differentiation with newton's method to minimize
#subdiagonals. Jacobian tends to become singular pretty quickly with random
#matrices
function [numDef] = msqriAD(H,s)
  MAXITER = 150;
  numDef = [];
  n = length(H);
  n1 = n;
  tol = sqrt(eps);
  shifts = eig(H(n-s+1:n,n-s+1:n));
  shifts = 100*randn(s,1);
  seed = eye(s);
  fprintf('starting msqriAD loop\n');
  do
    
    n = length(H);
    numDef = [numDef, n1 - n];
    fprintf('%d ',length(numDef));
    [B, ~,~,~,INFO] = trainBust(H,real(shifts),seed,false,false,'b',false);
    
    #check for deflation
    for i=n-s+1:n
      if(abs(B(i,i-1))/tol < abs(B(i,i)) + abs(B(i-1,i-1)))#can deflate
        H = H(1:i-1,1:i-1);
        for j = 1:s
          INFO.dots{j}.H = INFO.dots{j}.H(1:i-1,1:i-1);
        end
        break;
      end
    end
  
    #extract f(s) and J
    f = diag(H,-1);
    f = f(length(f) - s + 1:end);
    Jac = zeros(s);
    for i=1:s
      tmp = diag(INFO.dots{i}.H,-1);
      Jac(:,i) = tmp(length(tmp) - s + 1:end);
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

  until(length(H) < n1/2 || length(numDef) > MAXITER)
end
