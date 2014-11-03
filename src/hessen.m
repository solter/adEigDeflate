#[H] = hessen(A,toplt=false,top='b')
#forms hessenberg matrix. If hermitian, forms tridiagonal system.
#
function [H] = hessen(A, toplt = false, top = 'b')

  H = A;
  sym = false;
  if(norm(A - A') < sqrt(eps))
    sym = true;
  end
  %move bulge 1-column over per i
  for i=1:length(H)-2
    n = length(H);%min(i + bs,length(H));
    if(top == 'b')
      %create house vector
      v = H(i+1:n, i);
      v(1) += sgn(v(1))*sqrt(v'*v);
      v /= sqrt(v'*v);%normalize house vector
      
      %apply householder transformation to the right bits
      H(i+1:n,:) -= v*((2*v')*H(i+1:n,:));
      H(1:min(n+1,length(H)),i+1:n) -= (H(1:min(n+1,length(H)),i+1:n)*(2*v))*v';
      
      %zero out appropriate bits
      H(i+2:n,i) = 0;%clean up to ensure 0
      if(sym)
        H(i,i+2:n) = 0;%clean up to ensure 0
      end
    else
      %create house vector
      v = H(n-i+1,1:n-i);
      v(end) += sgn(v(end))*sqrt(v*v');
      v /= sqrt(v*v');
      
      %apply householder transformation to the right bits
      H(:,1:n-i) -= H(:,1:n-i)*(2*v')*v;
      H(1:n-i,:) -= v'*((2*v)*H(1:n-i,:));
      
      %zero out appropriate bits
      H(n-i+1,1:n-i-1) = 0;%clean up to ensure 0
      if(sym)
        H(1:n-i-1,n-i+1) = 0;%clean up to ensure 0
      end;
    end%if
    
    if(toplt)
      pltMat(H);
      pause(0.0);
    end
  end%for
end

function [sn] = sgn(v)
  if(v >= 0)
    sn= 1;
  else
    sn= -1;
  endif
endfunction

function [] = pltMat(H)
    imagesc(logNAN10(abs(H)));
    daspect([1 1 1]);
    h = colorbar;
    ytick = get(h, "ytick");
    set (h, "yticklabel", sprintf ('10^{%g}|', ytick));
end#function

function [out] = logNAN10(in)
  out = log10(in);
  out(isinf(out)) = nan;
  out(out < -16) = nan;
endfunction
