function [Q,R,ap] = swapSchur(Q,R,middle,_RELTOL);

% SYNTAX: [Q,R,ap] = SRSchur(Q,R,middle,_RELTOL)
%
% INPUT: orthogonal real Q and quasi-triangular real R such that AQ=QR.
% middle is the character 't', 'b', or 'm' for where the train was driven.
% _RELTOL is the tolerance for which deflation can occur 
%
% OUTPUT: orthogonal real Q and quasi-triangular real R such that AQ=QR 
% with good spikes
% 
% SUBFUNCTIONS: normalize.m, swaplist.m, select.m, swap.m, lu_complpiv.m
%
% SEE ALSO: schur.m, rsf2csf.m

r = find(abs(diag(R,-1)) > 100*eps);
s = 1:size(R,1)+1;
s(r+1) = [];
evals = length(s) - 1;

ap = [];
tGood = 0;%The rightmost tGood provide good spikes
bGood = 0;%The leftmost bGood provide good spikes
noGood = 0;%The number of evals which don't work anywhere

%save all the good spike values
while(checkT(Q,s,evals - tGood,middle) < _RELTOL); tGood++; end
while(checkB(Q,s,bGood+1,middle) < _RELTOL); bGood++; end

bQ = Q;
bR = R;
bs = s;

mIdx = ceil(length(Q)/2);
while(tGood + bGood + noGood < evals)
  for k = bGood + 1:evals - tGood -1
    %swap k'th eval with (k+1)st
    v      = s(k):s(k+1)-1; 
    w      = s(k+1):s(k+2)-1;
    nrA    = norm(R([v,w],[v,w]),inf);
    [Q,R]  = swap(Q,R,v,w); 
    s(k+1) = s(k)+s(k+2)-s(k+1);
    v      = s(k):s(k+1)-1; 
    w      = s(k+1):s(k+2)-1;
    if length(v)==2
      [Q,R] = normalize(Q,R,v);
    end
    if length(w)==2
      [Q,R] = normalize(Q,R,w);
    end
    ap  = [ap, norm(R(w,v),inf)/(10*eps*nrA)];

    %check to see if this is the best Q and R yet found, only considering
    %the outer half of the tip, weighted to push smaller values to outside
    if(middle == 'm')
      if( abs(Q(end,1:mIdx))*(mIdx:-1:1)' + abs(Q(1,mIdx:end))*(1:length(Q)-mIdx+1)' < ...
         abs(bQ(end,1:mIdx))*(mIdx:-1:1)' + abs(bQ(1,mIdx:end))*(1:length(Q)-mIdx+1)' 
      )
        bQ = Q; bR = R; bs = s;
      end
    elseif(middle == 't')
      %small values are near tips
      if( abs(Q(end,1:mIdx))*(mIdx:-1:1)' < abs(bQ(end,1:mIdx))*(mIdx:-1:1)' )
        bQ = Q; bR = R; bs = s;
      end
    else
      if( abs(Q(1,mIdx:end))*(1:length(Q)-mIdx+1)' < abs(bQ(1,mIdx:end))*(1:length(Q)-mIdx+1)' )
        bQ = Q; bR = R; bs = s;
      end
    end

  end
  bg = checkB(Q,s,bGood+1,middle);
  tg = checkT(Q,s,evals - tGood,middle);
  bGood += bg < _RELTOL;
  tGood += tg < _RELTOL;
  noGood += tg > _RELTOL;%if the one pushed to the end didn't work
end

%if not deflatable, use best one found
if((middle ~= 'm' && bGood +tGood == 0) ||...
    bGood + tGood < mIdx)
  Q = bQ;
  R = bR;
  s = bs;
end

R = R - tril(R,-2);
for j=2:length(s)-1; R(s(j),s(j)-1)=0; end

%--------------------------------------------%
%functions to check the spikes for zeros at top and bottom
%return true if there's a zero there
  function [b] = checkT(Q,s,tIdx,middle)
    if(middle == 't')%if we don't use this row, don't check it
      b = inf;
    else
      %this is the size of the entry for this eval.
      %if complex conjugate eval, then takes the max.
      b = max(abs(Q(1,s(tIdx):s(tIdx+1)-1)));
    end

  function [b] = checkB(Q,s,bIdx,middle)
    if(middle == 'b')
      b = inf;
    else
      b = max(abs(Q(end,s(bIdx):s(bIdx+1)-1)));
    end

% ----------------------------------------------%

  function [U,S] = normalize(U,S,v);
  %normalizes complex conjugate blocks
  n  = size(S,1);
  Q  = rot(S(v,v));%givens rotation
  S(:,v) = S(:,v)*Q;
  S(v,:) = Q'*S(v,:); 
  U(:,v) = U(:,v)*Q;

  % ----------------------------------------------%

  function Q = rot(X);
  c = 1; s = 0;
  if X(1,1)~=X(2,2);
    tau   = (X(1,2)+X(2,1))/(X(1,1)-X(2,2));
    off   = sqrt(tau^2+1);
    v     = [tau - off, tau + off];
    [d,w] = min(abs(v));
    c     = 1/sqrt(1+v(w)^2);
    s     = v(w)*c;
  end
  Q = [c -s;s c];

% -----------------------------------------------%

  function [U,S] = swap(U,S,v,w);  
  [p,q] = size(S(v,w)); Ip = eye(p); Iq = eye(q);
  r = [];
  for j=1:q
    r = [r;S(v,w(j))];
  end 
  K = kron(Iq,S(v,v))-kron(S(w,w)',Ip);
  [L,H,P,Q] = lu_complpiv(K);
  e = min(abs(diag(H)));
  sigp = 1:p*q;
  for k = 1:p*q-1;
    sigp([k,P(k)]) = sigp([P(k),k]);
  end 
  r = e*r(sigp); 
  x = (H\(L\r));
  sigq = 1:p*q;
  for k = 1:p*q-1;
    sigq([k,Q(k)]) = sigq([Q(k),k]);
  end
  x(sigq) = x;
  X = [];
  for j=1:q
    X = [X,x((j-1)*p+1:j*p)];
  end
  [Q,R]      = qr([-X;e*Iq]);
  S(:,[v,w]) = S(:,[v,w])*Q;
  S([v,w],:) = Q'*S([v,w],:); 
  U(:,[v,w]) = U(:,[v,w])*Q;

  % -----------------------------------------------%

    function [L,U,P,Q] = lu_complpiv(A);
    P = []; Q = []; n = size(A,1);
    for k=1:n-1;
      [a,r] = max(abs(A(k:n,k:n)));
      [dummy,c] = max(abs(a));
      cl  = c+k-1;
      rw  = r(c)+k-1;
      A([k,rw],:) = A([rw,k],:);
      A(:,[k,cl]) = A(:,[cl,k]);
      P(k) = rw; Q(k) = cl;
      if A(k,k) ~= 0;
        rs = k+1:n;
        A(rs,k)  = A(rs,k)/A(k,k);
        A(rs,rs) = A(rs,rs)-A(rs,k)*A(k,rs);
      end
    end
    U = tril(A')'; L = tril(A,-1) + eye(n);
   
  

  
