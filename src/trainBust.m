##-*- texINFO -*-
##@deftypefn{Function File} {[@var{H}, @var{INFO}] =} trainBust(@var{A},@var{shifts})
##@deftypefnx{Function File} {[@var{H}, @var{INFO}] =} trainBust(@var{A},@var{shifts},@var{toplt})
##@deftypefnx{Function File} {[@var{H}, @var{INFO} =} trainBust(@var{A},@var{shifts},@var{toplt},@var{toprt})
##@deftypefnx{Function File} {[@var{H}, @var{INFO} =} trainBust(@var{A},@var{shifts},@var{toplt},@var{toprt},@var{middle})
##
##This introduces 2 trains of bulges (1 from top and 1 from bottom) which sweep until they meet.@*
##Then a Schur decomposition is performed to create spikes replacing the bulges.@*
##
##Inputs:@*
##  @var{A} - A hessenberg matrix. Results will be meaningless if not hessenberg.@*
##  @var{shifts} - Lists of shifts to introduce in the bulges.
##    For indentical shifts in top and bottom, this should be symmetric.
##    AKA: @var{shifts}(1) = @var{shifts}(end); @var{shifts}(2) = @var{shifts}(end - 1); ...@*
##  @var{shiftSeed} - Matrix of directions for shift derivatives.
##    Each row should represent a different direction, thus for d directions and s
##    shifts this matrix will be d by s.@*
##  @var{toplt} - Optional argument. If it is true, this function
##    will plot each step of the bulge creation and chasing.@*
##  @var{toprt} - Optional argument. If it is true, this function
##    will create a directory ../impStepPlts and save the plots to it.@*
##  @var{middle} - If 'm' drives two train bulges to center and spike it (default).
##    If 'b' will drive a train bulge to bottom and spike it.
##    If 't' will drive a train bulge to top and spike it.
##
##Outputs:@*
##  @var{H} - The matrix A after running bulge trains into it's center and
##    creating the spikes@*
##  @var{shiftDerivs} - The derivative of the spikes w.r.t. the shifts in the directions
##    specified by @var{shiftSeed}. This will only represent the smallest values
##    in the spike, so will be the same dimensions as @var{shiftSeed}.
##  @var{INFO} - a structure with INFOrmation. Fields:
##    vspIdx -> vertical spike's column index
##    hspIdx -> horizontal spike's row index (only for middle)
##    nz -> count of zeros in spikes
## @end deftypefn
function [H, shiftDerivs, INFO] = trainBust(H, Shifts, shiftSeed, toplt = false, toprt = false, middle = 'm')

  
  A = H;
  #TODO: actually deflate matrix, so far only finds deflation points while pushing bulges
  #need to do deflation logic for spike

  #TODO 2: implement derivatives... use in deflation logic

  _PAUSELEN = 0.0;#pause time when playing real time
  _MAXBULGESIZE = 5;#maximum bulge size
  _EXTRASPIKE = 0;#floor(.5*length(Shifts));#extra size of spike beyond shifts
  _RELTOL = 1e-6;
  n = length(H);#convenience shorthand
  toprt = toprt && toplt;%must be plotting to print
  INFO = struct();

  if(toprt)
    mkdir('../impStepPlts');
    pltNum = 0;
  else
    toprt = false;
  end#if

  tstIdx = [];#top bulge starts
  tendIdx = -1;
  bstIdx = [];#bottom bulge starts
  bendIdx = -1;
  do
    #TOP BULGE LOGIC
    ++tstIdx;
    ++tendIdx;

    if(middle != 't')
      #move top bulges over
      tempEndIdx = tendIdx;
      for tempStIdx = tstIdx#for each
        %perform householder step
        H = hhStepR(H,tempStIdx,tempEndIdx+1);
        tempEndIdx = tempStIdx;
      endfor

      #add shift to top if still need to add shifts
      if(middle == 'b')
        #The +1 is because the bulge will contain bsize - 1 shifts
        bsize = min( _MAXBULGESIZE,length(Shifts)) + 1;
      elseif(middle == 'm')
        #The +1 is because the bulge will contain bsize - 1 shifts
        bsize = min( _MAXBULGESIZE, ceil(length(Shifts)/2)) + 1;
      end
      if(!isempty(Shifts) && (isempty(tstIdx) || tstIdx(end) > bsize ) );
        #calculate necessary polynomial bits
        polyH = eye(bsize);
        #explicitly form the polynomial of shifts to shove in.
        #Note that it should be 1:bsize-1 to account for fill in
        for i=1:bsize-1
          polyH *= H(1:bsize,1:bsize) - Shifts(1)*eye(bsize);
          Shifts(1) = [];
        endfor

        #create house vector
        H = hhStepR(H,1,bsize,polyH(:,1));
        tstIdx = [tstIdx 1];#add this to index (shift) table
        if(tendIdx == 0)
          tendIdx = bsize;
        end
      endif
    endif

    #BOTTOM TRAIN
    if(middle == 't' || (middle == 'm' && tendIdx + bendIdx + 1 < n))#if they won't intersect
      ++bendIdx;
      ++bstIdx;
       
      #shift up the bottom bulges
      tempEndIdx = bendIdx;
      for tempStIdx = bstIdx#for each
        H = hhStepB(H,tempStIdx,tempEndIdx);
        tempEndIdx = tempStIdx;
      endfor
      
      #add shift to bottom if still need to add shifts
      #The +1 is because the bulge will contain bsize - 1 shifts
      bsize = min( _MAXBULGESIZE,length(Shifts))+1;
      if(!isempty(Shifts) && (isempty(bstIdx) || bstIdx(end) + 1 > bsize ) )
        #calculate necessary polynomial bits
        polyH = eye(bsize);
        #explicitly form the polynomial of shifts to shove in.
        #Note that it should be 1:bsize-1 to account for fill in
        for i=1:bsize-1
          polyH *= H(n-bsize+1:n,n-bsize+1:n) - Shifts(end)*eye(bsize);
          Shifts(end) = [];
        endfor
        
        H = hhStepB(H,0,bsize-1,polyH(end,:));
        #create house vector
        bstIdx = [bstIdx 0];#add this to shift table
        if(bendIdx == 0)
          bendIdx = bsize;
        end
      endif
      
    endif

    if(toplt)
      pltMat(H);
      if(toprt)
        print(sprintf('../impStepPlts/impstep%03d.png',++pltNum));
      else
        pause(_PAUSELEN);
      endif
    end#if
    fflush(stdout);

    doneb = middle == 'b' && tendIdx + 1 >= n;%check for hitting bottom
    donet = middle == 't' && bendIdx + 1 >= n;%check for hitting bottom
    donem = middle == 'm' && tendIdx + bendIdx + 3 > n;%check for meeting at middle
  until(donet || donem || doneb)#if the trains are touching

  INFO.nz = 0;

  if(middle == 'm')#if doing middle deflation
    spSt = tstIdx(end);
    spEnd = n - bstIdx(end);
  elseif(middle == 'b')#if doing bottom deflation 
    spSt = tstIdx(end) - _EXTRASPIKE;#index of block
    spEnd = n;
  elseif(middle == 't')#if doing top deflation 
    spSt = 1;
    spEnd = n - bstIdx(end) + 1 + _EXTRASPIKE;#index of block
  endif
    
  [H, INFO] = schurComp(H, spSt, spEnd, INFO, _RELTOL);

  if(toplt)
    pltMat(H);
    if(toprt)
      print(sprintf('../impStepPlts/impstep%03d.png',++pltNum));
    else
      pause(_PAUSELEN);
    endif
  endif

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

#perform a householder transformation, pushing a bulge rightwards
function [H] = hhStepR(H,stIdx, eIdx,v=[])

  #create house vector
  if(length(v) ~= eIdx - stIdx + 1)#if a vector is not provided
    v = H(stIdx:eIdx, stIdx-1);
  end
  v(1) += sgn(v(1))*sqrt(v'*v);
  v /= sqrt(v'*v);#normalize house vector

  #apply householder transformation to the right bits
  H(stIdx:eIdx,:) -= v*((2*v')*H(stIdx:eIdx,:));
  maxEff = min(length(H),eIdx + 1);
  H(1:maxEff,stIdx:eIdx) -= (H(1:maxEff,stIdx:eIdx)*(2*v))*v';
  if(stIdx > 1)
    H(stIdx+1:maxEff,stIdx-1) = 0;#zeros everything out for exactness
  end
end

#perform a householder transformation, pushing a bulge from bottom
function [H] = hhStepB(H,stIdx, eIdx,v=[])
  n = length(H);
  
  #create house vector
  if(length(v) ~= eIdx - stIdx + 1)#if a vector is not provided
    v = H(n - stIdx + 1, n - eIdx:n - stIdx);
  end
  v(end) += sgn(v(end))*sqrt(v*v');
  v /= sqrt(v*v');#normalize house vector
  
  #apply householder transformation to the right bits
  minEff = max(1,n - eIdx - 1);
  H(n - eIdx:n - stIdx, minEff:n) -= v'*((2*v)*H(n - eIdx:n - stIdx,minEff:n));
  H(:,n - eIdx:n - stIdx) -= (H(:,n - eIdx:n - stIdx)*(2*v'))*v;
  if(stIdx > 0)
    H(n - stIdx + 1,minEff: n - stIdx-1) = 0;#zeros everything out for exactness
  end
end

#perform the schur decomposition
function [H,INFO] = schurComp(H,spSt,spEnd,INFO,_RELTOL)
  
  n = length(H);
  INFO.nz = 0;

  #actually perform the schur decomposiotn
  [spRot, H(spSt:spEnd,spSt:spEnd)] = schur(H(spSt:spEnd,spSt:spEnd));
  if(spSt > 1)%apply it to the left/top portion
    H(1:spSt-1,spSt:spEnd) = H(1:spSt-1,spSt:spEnd) * spRot;
  end
  if(spEnd < n)%apply it to the right/bottom portion
    H(spSt:spEnd,spEnd+1:n) = spRot'*H(spSt:spEnd,spEnd+1:n);
  end

  #form spikes, make hard zeros, and put into H
  vSpike = spRot'(:,1);
  hSpike = spRot(end,:);
  vsz = abs(vSpike) < _RELTOL;
  hsz = abs(hSpike) < _RELTOL;
  vSpike(vsz) = 0;
  hSpike(hsz) = 0;
  
  
  if(spEnd < n)%apply it to the right/bottom portion
    H(spEnd+1,spSt:spEnd) = H(spEnd+1,spEnd) * hSpike;
    INFO.nz += sum(hsz);
    INFO.hspIdx = spEnd + 1;
  end
  if(spSt > 1)%apply it to the left/top portion
    H(spSt:spEnd,spSt-1) = H(spSt,spSt-1) * vSpike;
    INFO.nz += sum(vsz);
    INFO.vspIdx = spSt - 1;
  end

end
