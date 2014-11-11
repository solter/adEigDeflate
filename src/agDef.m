function [numDef] = agDef(H,s)
  MAXITER = 150;
  numDef = [];
  n = length(H);
  n1 = n;
  tol = sqrt(eps);
  do
    n = length(H);
    numDef = [numDef, n1 - n];
    shifts = eig(H(n-s+1:n,n-s+1:n));
    H = trainBust(H,real(shifts),[],false,false,'b',true);
    n = length(H);
    
    #check for deflation
    for i=n-s+1:n
      if(abs(H(i,i-1))/tol < abs(H(i,i)) + abs(H(i-1,i-1)))#can deflate
        H = H(1:i-1,1:i-1);
        break;
      end
    end
    
  until(length(H) < n1/2 || length(numDef) > MAXITER )
end
