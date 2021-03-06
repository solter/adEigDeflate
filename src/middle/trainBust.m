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
##  @var{shifts}, - Lists of shifts to introduce in the bulges.
##    For indentical shifts in top and bottom, this should be symmetric.
##    AKA: @var{shifts}(1) = @var{shifts}(end); @var{shifts}(2) = @var{shifts}(end - 1); ...@*
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
##  @var{INFO} - a structure with INFOrmation. Fields:
##    vspIdx -> vertical spike's column index
##    hspIdx -> horizontal spike's row index (only for middle)
##    nz -> count of zeros in spikes
## @end deftypefn
function [H, INFO] = trainBust(H, Shifts, toplt = false, toprt = false, middle = 'm')
  A = H;
  #TODO: actually deflate matrix, so far only finds deflation points while pushing bulges
  #need to do deflation logic for spike

  #TODO 2: implement derivatives... use in deflation logic

  _PAUSELEN = 0.0;#pause time when playing real time
  _MAXBULGESIZE = 4;#maximum bulge size
  _EXTRASPIKE = 0;#floor(.5*length(Shifts));#extra size of spike beyond shifts
  _RELTOL = 1e-6;
  n = length(H);#convenience shorthand
  toprt = toprt && toplt;%must be plotting to print

  if(toprt)
    mkdir('../impStepPlts');
    pltNum = 0;
  else
    toprt = false;
  end#if

  tstIdx = [];#top bulge starts
  bstIdx = [];#bottom bulge starts
  #end of bottom bulge on top train
  tendIdx = min(_MAXBULGESIZE,ceil(length(Shifts)/2));
  #how far top of bottom train is from the right side
  bendIdx = min(_MAXBULGESIZE,floor(length(Shifts)/2)) - 1;
  do
    #TOP BULGE LOGIC
    ++tstIdx;
    ++tendIdx;

    if(middle != 't')
      #move top bulges over
      tempEndIdx = tendIdx;
      for tempStIdx = tstIdx#for each
        #create house vector
        v = H(tempStIdx:tempEndIdx, tempStIdx-1);
        v(1) += sgn(v(1))*sqrt(v'*v);
        v /= sqrt(v'*v);#normalize house vector
        #apply householder transformation to the right bits
        H(tempStIdx:tempEndIdx,:) -= v*((2*v')*H(tempStIdx:tempEndIdx,:));
        H(1:tempEndIdx+1,tempStIdx:tempEndIdx) -= (H(1:tempEndIdx+1,tempStIdx:tempEndIdx)*(2*v))*v';
        H(tempStIdx+1:tempEndIdx,tempStIdx-1) = 0;#zeros everything out for exactness
        tempEndIdx = tempStIdx;
      endfor

      #add shift to top if still need to add shifts
      bsize = min( _MAXBULGESIZE,ceil(length(Shifts)/2))+1;
      if(!isempty(Shifts) && (isempty(tstIdx) || tstIdx(end) > bsize ) )
        #calculate necessary polynomial bits
        polyH = eye(bsize);
        for i=1:bsize-1
          polyH *= H(1:bsize,1:bsize) - Shifts(1)*eye(bsize);
          Shifts(1) = [];
        endfor
        
        #create house vector
        v = polyH(1:bsize, 1);
        v(1) += sgn(v(1))*sqrt(v'*v);
        v /= sqrt(v'*v);#normalize house vector
        #apply householder transformation to the right bits
        H(1:bsize,:) -= v*((2*v')*H(1:bsize,:));
        H(1:bsize+1,1:bsize) -= (H(1:bsize+1,1:bsize)*(2*v))*v';#many zero mulitplies
        tstIdx = [tstIdx 1];#add this to shift table
      endif
    endif

    #BOTTOM TRAIN
    ++bendIdx;
    ++bstIdx;
    if(middle == 't' || (middle == 'm' && tendIdx + bendIdx + 1 < n))#if they won't intersect
    
      #shift up the bottom bulges
      tempEndIdx = bendIdx;
      for tempStIdx = bstIdx#for each
        #create house vector
        v = H(n - tempStIdx + 1, n - tempEndIdx:n - tempStIdx);
        v(end) += sgn(v(end))*sqrt(v*v');
        v /= sqrt(v*v');#normalize house vector
        #apply householder transformation to the right bits
        H(n - tempEndIdx:n - tempStIdx,n - tempEndIdx - 1:n) -= v'*((2*v)*H(n - tempEndIdx:n - tempStIdx,n - tempEndIdx - 1:n));
        H(:,n - tempEndIdx:n - tempStIdx) -= (H(:,n - tempEndIdx:n - tempStIdx)*(2*v'))*v;
        H(n - tempStIdx + 1,n - tempEndIdx: n - tempStIdx-1) = 0;#zeros everything out for exactness
        tempEndIdx = tempStIdx;
      endfor

      #add shift to bottom if still need to add shifts
      bsize = min( _MAXBULGESIZE,length(Shifts))+1;
      if(!isempty(Shifts) && (isempty(bstIdx) || bstIdx(end) + 1 > bsize ) )
        #calculate necessary polynomial bits
        polyH = eye(bsize);
        for i=1:bsize-1
          polyH *= H(n-bsize+1:n,n-bsize+1:n) - Shifts(end)*eye(bsize);
          Shifts(end) = [];
        endfor
        
        #create house vector
        v = polyH(end,:);
        v(end) += sgn(v(end))*sqrt(v*v');
        v /= sqrt(v*v');#normalize house vector
        #apply householder transformation to the right bits
        H(n-bsize+1:n,n - bsize - 1:n) -= v'*((2*v)*H(n-bsize+1:n,n - bsize - 1:n));
        H(:,n-bsize+1:n) -= (H(:,n-bsize+1:n)*(2*v'))*v;#many zero mulitplies
        bstIdx = [bstIdx 0];#add this to shift table
      endif

    else
      #reset these since nothing was done
      --bendIdx;
      --bstIdx;
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
    donet = middle == 't' && bendIdx + 2 >= n;%check for hitting bottom
    donem = middle == 'm' && tendIdx + bendIdx + 3 > n;%check for meeting at middle
  until(donet || donem || doneb)#if the trains are touching

  if(middle == 'm')#if doing middle deflation
    spSt = tstIdx(end);
    spEnd = n - bstIdx(end);
    #compute Schur decomp and apply it to right bits
    [spRot, H(spSt:spEnd,spSt:spEnd)] = schur(H(spSt:spEnd,spSt:spEnd));
    H(1:spSt-1,spSt:spEnd) = H(1:spSt-1,spSt:spEnd) * spRot;
    H(spSt:spEnd,spEnd+1:end) = spRot'*H(spSt:spEnd,spEnd+1:end);
    #form spikes, make hard zeros, and put into H
    vSpike = spRot'(:,1);
    hSpike = spRot(end,:);

    vsz = abs(vSpike) < _RELTOL;
    hsz = abs(hSpike) < _RELTOL;
    INFO.nz = sum(vsz) + sum(hsz);
    vSpike(vsz) = 0;
    hSpike(hsz) = 0;
    INFO.vspIdx = spSt - 1;
    INFO.hspIdx = spEnd+1;
    H(spSt:spEnd,spSt-1) = H(spSt,spSt-1) * vSpike;
    H(spEnd+1,spSt:spEnd) = H(spEnd+1,spEnd) * hSpike;
  elseif(middle == 'b')#if doing bottom deflation 
    spSt = tstIdx(end) - _EXTRASPIKE;#index of block
    [spRot,H(spSt:end,spSt:end)] = schur(H(spSt:end,spSt:end));%compute schur decomp and do lower right portion
    H(1:spSt-1,spSt:end) = H(1:spSt-1,spSt:end)*spRot;%upper right portion of H
    #spike and zero out crappy stuff
    sp = spRot'(:,1);
    spz = abs(sp) < _RELTOL;
    sp(spz) = 0;
    INFO.nz = sum(spz);
    INFO.vspIdx = spSt - 1;
    H(spSt:end,spSt-1) = H(spSt,spSt-1)*sp;
  elseif(middle == 't')#if doing top deflation 
    spEnd = n - bstIdx(end) + 1 + _EXTRASPIKE;#index of block
    [spRot,H(1:spEnd,1:spEnd)] = schur(H(1:spEnd,1:spEnd));%compute schur decomp and do upper left portion
    H(1:spEnd,spEnd+1:end) = spRot'*H(1:spEnd,spEnd+1:end);%upper right portion of H
    #spike and zero out crappy stuff
    sp = spRot(end,:);
    spz = abs(sp) < _RELTOL;
    sp(spz) = 0;
    INFO.nz = sum(spz);
    INFO.hspIdx = spEnd+1;
    H(spEnd + 1,1:spEnd) = H(spEnd + 1,spEnd)* sp;
  endif

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
