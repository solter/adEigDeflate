##-*- texinfo -*-
##@deftypefn{Function File} {[@var{H}, @var{spSt}, @var{spEnd}, @var{pltNum}] =} agErlyDef(@var{A},@var{topshifts})
##@deftypefnx{Function File} {[@var{H}, @var{spSt}, @var{spEnd}, @var{pltNum}] =} agErlyDef(@var{A},@var{topshifts},@var{toplt})
##
##This introduces 2 bulges (1 from top and 1 from bottom) which sweep until they meet.@*
##Then a Schur decomposition is performed to create spikes replacing the bulges.@*
##
##Inputs:@*
##  @var{A} - A hessenberg matrix. Results will be meaningless if not hessenberg.@*
##  @var{topshifts}, - Lists of shifts to introduce in the bulges respectively.
##    Should be less than ~6 to avoid shift blurring.@*
##  @var{toplt} - Optional argument. If it is true, this function
##    will plot each step of the bulge creation and chasing.@*
##  @var{toprt} - Optional argument. If it is true, this function
##    will create a directory ../impStepPlts save the plots to it.@*
##
##Outputs:@*
##  @var{H} - The matrix A after running 2 bulges into it's center and
##    creating the spikes@*
## @end deftypefn
function [H, pltNum] = bulgeChase(A, Shifts, toplt = false, toprt = false)
  H=A;

  pl = .05;#pause time when playing real time
  
  if(toplt)
    mkdir('../impStepPlts');
    pltNum = 0;
  else
    toprt = false;
  end#if

  stIdx = [];#bulge starts
  endIdx = 1;#end of last bulge
  bsize = 0;#bulge size
  do
    ++stIdx;
    ++endIdx;
    
    #move bulges over
    tendIdx = endIdx;
    for tstIdx = stIdx#for each 
      #create house vector
      v = H(tstIdx:tendIdx, tstIdx-1);
      v(1) += sgn(v(1))*sqrt(v'*v);
      v /= sqrt(v'*v);#normalize house vector
      #apply householder transformation to the right bits
      H(tstIdx:tendIdx,:) -= v*((2*v')*H(tstIdx:tendIdx,:));
      H(1:tendIdx+1,tstIdx:tendIdx) -= (H(1:tendIdx+1,tstIdx:tendIdx)*(2*v))*v';
      H(tstIdx+1:tendIdx,tstIdx-1) = 0;#zeros everything out for exactness
      tendIdx = tstIdx;
    endfor   

    #add shift to top if still need to add shifts
    if(!isempty(Shifts) && (isempty(stIdx) || stIdx(end) > 2))
      #create house vector
      v = H(1:tendIdx, 1);
      v(1) -= Shifts(1);
      Shifts(1) = [];
      v(1) += sgn(v(1))*sqrt(v'*v);
      v /= sqrt(v'*v);#normalize house vector
      #apply householder transformation to the right bits
      H(1:tendIdx,:) -= v*((2*v')*H(1:tendIdx,:));
      H(:,1:tendIdx) -= (H(:,1:tendIdx)*(2*v))*v';#many zero mulitplies
     
      if(++bsize > 3 || isempty(Shifts))
        stIdx = [stIdx 1];
        bsize = 0;
      endif
    endif

    if(toplt)
      pltMat(H);
      if(toprt)
        print(sprintf('../impStepPlts/impstep%03d.png',++pltNum));
      else
        pause(pl);
      endif
    end#if
    fflush(stdout);
  until(endIdx + 1 >= length(H))#if it touches the bottom 

  ++stIdx;

  do
    
    endIdx = length(H);
    for tstIdx = stIdx#for each 
      #create house vector
      v = H(tstIdx:endIdx, tstIdx-1);
      v(1) += sgn(v(1))*sqrt(v'*v);
      v /= sqrt(v'*v);#normalize house vector
      #apply householder transformation to the right bits
      H(tstIdx:endIdx,:) -= v*((2*v')*H(tstIdx:endIdx,:));
      tendIdx = min(endIdx + 1, length(H));
      H(1:tendIdx,tstIdx:endIdx) -= (H(1:tendIdx,tstIdx:endIdx)*(2*v))*v';
      H(tstIdx+1:endIdx,tstIdx-1) = 0;#zeros everything out for exactness
      endIdx = tstIdx;
    endfor
    
    ++stIdx;
    if(stIdx(1) >= length(H))
      stIdx(1) = [];
    endif
    
    if(toplt)
      pltMat(H);
      if(toprt)
        print(sprintf('../impStepPlts/impstep%03d.png',++pltNum));
      else
        pause(pl);
      endif
    end#if

  until(isempty(stIdx))

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
