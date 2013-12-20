%H must be hessenberg, shifts should be a vector
%OUTPUTS: H - 1-sweep after H
function [H] = hessen(A)

  H = A;
  %move bulge 1-column over per i
  for i=1:length(H)-2
    endIdx = length(H);%min(i + bs,length(H));
    %create house vector
    v = H(i+1:endIdx, i);
    if(v(1) > 0)
      v(1) += sqrt(v'*v);
    else
      v(1) -= sqrt(v'*v);
    endif
    v /= sqrt(v'*v);%normalize house vector
    
    %apply householder transformation to the right bits
    H(i+1:endIdx,:) -= v*((2*v')*H(i+1:endIdx,:));
    H(1:min(endIdx+1,length(H)),i+1:endIdx) -= (H(1:min(endIdx+1,length(H)),i+1:endIdx)*(2*v))*v';
    %zero out appropriate bits
    H(i+2:length(H),i) = 0;%clean up to ensure 0
end%function
