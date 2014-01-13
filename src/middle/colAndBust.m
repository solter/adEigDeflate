##-*- texinfo -*-
##@deftypefn{Function File} {[@var{H}, @var{spSt}, @var{spEnd}, @var{pltNum}] =} colAndBust(@var{A},@var{topshifts}, @var{botShifts})
##@deftypefnx{Function File} {[@var{H}, @var{spSt}, @var{spEnd}, @var{pltNum}] =} colAndBust(@var{A},@var{topshifts}, @var{botShifts},@var{toplt})
##
##This introduces 2 bulges (1 from top and 1 from bottom) which sweep until they meet.@*
##Then a Schur decomposition is performed to create spikes replacing the bulges.@*
##
##Inputs:@*
##  @var{A} - A hessenberg matrix. Results will be meaningless if not hessenberg.@*
##  @var{topshifts}, @var{botshifts} - Lists of shifts to introduce in the top and bottom bulges respectively.
##    Should be less than ~6 to avoid shift blurring.@*
##  @var{toplt} - Optional argument. If it is true, this function
##    will create a directory ../impStepPlts and plot each step of the
##    bulge creation and chasing to it.@*
##
##Outputs:@*
##  @var{H} - The matrix A after running 2 bulges into it's center and
##    creating the spikes@*
##  @var{spSt}, @var{spEnd} - the indices the spike's column and row respectively
## @end deftypefn
function [H, spSt, spEnd, pltNum] = colAndBust(A, topShifts, botShifts, toplt=false, toprt=false)
  H=A;

  tsl = length(topShifts);
  bsl = length(botShifts);
  _PAUSELEN = 0.1;#length to pause when printing

    if(toplt)
      pltMat(H);
      if(toprt)
        pltNum=0;
        print(sprintf('../impStepPlts/impstep%03d.png',++pltNum));
      else
        pause(_PAUSELEN);
      endif
    end#if
  #put shifts into H

  #add bulge to the top
  for i=1:tsl
    endIdx = 1+i;

    #create house vector
    v = H(1:endIdx, 1);
    v(1) -= topShifts(i);
    v(1) += sgn(v(1))*sqrt(v'*v);
    v /= sqrt(v'*v);#normalize house vector
    #apply householder transformation to the right bits
    H(1:endIdx,:) -= v*((2*v')*H(1:endIdx,:));
    H(:,1:endIdx) -= (H(:,1:endIdx)*(2*v))*v';#many zero mulitplies

    if(toplt)
      pltMat(H);
      if(toprt)
        print(sprintf('../impStepPlts/impstep%03d.png',++pltNum));
      else
        pause(_PAUSELEN);
      endif
    end#if
  endfor

  #add bulge to the bottom
  for i=1:bsl
    endIdx = 1+i;
    v = H(end,(end-endIdx):end);
    v(end) -= botShifts(i);
    v(end) += sgn(v(end))*sqrt(v*v');
    v /= sqrt(v*v');#normalize house vector
    H((end-endIdx):end,:) -= v'*((2*v)*H((end-endIdx):end,:));#many zero multiplies
    H(:,(end-endIdx):end) -= (H(:,(end-endIdx):end)*(2*v'))*v;

    if(toplt)
      pltMat(H);
      if(toprt)
        print(sprintf('../impStepPlts/impstep%03d.png',++pltNum));
      else
        pause(_PAUSELEN);
      endif
    end#if
  end#for

  #move bulge 1-column over per i
  n = length(H);
  tstIdx = 1;
  tendIdx = 1 + tsl;
  bendIdx = n;
  bstIdx = n - bsl;
  do
    ++tstIdx;
    ++tendIdx;

    #move top bulge
    #create house vector
    v = H(tstIdx:tendIdx, tstIdx-1);
    v(1) += sgn(v(1))*sqrt(v'*v);
    v /= sqrt(v'*v);#normalize house vector
    #apply householder transformation to the right bits
    H(tstIdx:tendIdx,:) -= v*((2*v')*H(tstIdx:tendIdx,:));
    H(1:tendIdx+1,tstIdx:tendIdx) -= (H(1:tendIdx+1,tstIdx:tendIdx)*(2*v))*v';
    H(tstIdx+1:tendIdx,tstIdx-1) = 0;#zeros everything out for exactness

    --bendIdx;
    --bstIdx;
    if(tendIdx + 1 < bstIdx)#if they won't intersect
      #move bottom bulge
      v = H(bendIdx + 1,bstIdx:bendIdx);
      v(end) += sgn(v(end))*sqrt(v*v');
      v /= sqrt(v*v');#normalize house vector
      H(bstIdx:bendIdx,bstIdx - 1:n) -= v'*((2*v)*H(bstIdx:bendIdx,bstIdx-1:n));
      H(:,bstIdx:bendIdx) -= (H(:,bstIdx:bendIdx)*(2*v'))*v;
      H(bendIdx + 1,bstIdx:bendIdx-1) = 0;#zeros stuff out for exactness
    else
      bendIdx++;
      bstIdx++;;
    end#if

    if(toplt)
      pltMat(H);
      if(toprt)
        print(sprintf('../impStepPlts/impstep%03d.png',++pltNum));
      else
        pause(_PAUSELEN);
      endif
    end#if
  until(tendIdx + 2 >= bstIdx)#if they are touching


  #create spikes
  spSt = tstIdx+1;#index of left block
  spEnd = bendIdx-1;
  [spRot,~] = schur(H(spSt:spEnd,spSt:spEnd));
  H(1:spEnd+1,spSt:spEnd) = H(1:spEnd+1,spSt:spEnd)*spRot;
  H(spSt:spEnd,(spSt-1):n) = spRot'*H(spSt:spEnd,(spSt-1):n);

  spSt--;
  spEnd++;

    if(toplt)
      pltMat(H);
      if(toprt)
        print(sprintf('../impStepPlts/impstep%03d.png',++pltNum));
      else
        pause(_PAUSELEN);
      endif
    end#if

  if(toprt)
    close all;
  end#if

end#function

function [sn] = sgn(v)
  if(v >= 0)
    sn= 1;
  else
    sn= -1;
  endif
endfunction

function [] = pltMat(H)
    imagesc(log(abs(H)) + .01);
    daspect([1 1 1]);
end#function
