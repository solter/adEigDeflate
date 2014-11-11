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
##    creating the spikes.
##    If matrix was split at middle, this will be a cell array with 2 matrices
#     which were split from the split matrix H.@*
##  @var{egs} - eigenvalues which have been deflated out of matrix
##  @var{shiftDerivs} - The derivative of the spikes w.r.t. the shifts in the directions
##    specified by @var{shiftSeed}. 
##  @var{INFO} - a structure with INFOrmation. Fields:
##    vspIdx -> vertical spike's column index
##    hspIdx -> horizontal spike's row index (only for middle)
##    nz -> count of zeros in spikes
## @end deftypefn
function [H, egs, shiftDerivs, spikes, INFO] = trainBust(H, Shifts, shiftSeed, toplt = false, toprt = false, middle = 'm', toSpike = true)

  _PAUSELEN = 0.0;#pause time when playing real time
  _MAXBULGESIZE = 5;#maximum bulge size
  _EXTRASPIKE = 0;#floor(.5*length(Shifts));#extra size of spike beyond shifts
  _RELTOL = sqrt(eps);
  n = length(H);#convenience shorthand
  toprt = toprt && toplt;#must be plotting to print
  INFO = struct();
  sym = norm(H - H','inf') < sqrt(eps);
  %Make sure the shifts are real.
  %Complex shifts broke this...
  %in that shoving them in the matrix was no longer similar to the initial matrix
  Shifts = real(Shifts);
  shiftSeed = real(shiftSeed);
  if(middle == 'm'); toSpike = true; end

  #the structure to hold all the derivative information
  #make sure shiftSeed is flat
  if(~isempty(shiftSeed) && size(shiftSeed,2) ~= length(Shifts))
    shiftSeed = shiftSeed';
    if(size(shiftSeed,2) ~= length(Shifts))
      error("one dimension of shiftSeed should be number of shifts");
    end
  end

  s = size(shiftSeed,1);
  dots = cell(s,1);

  for i=1:s
    dots{i}.seed = shiftSeed(i,:);
    dots{i}.H = zeros(size(H));
  end
  
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
        [H, dots] = hhStepR(H,tempStIdx,tempEndIdx+1,dots,sym);
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
        for j=1:s
          dots{j}.polyH = zeros(bsize);
        end
        polyH = eye(bsize);
        #explicitly form the polynomial of shifts to shove in.
        #Note that it should be 1:bsize-1 to account for fill in
        for i=1:bsize-1
          for j=1:s
          dots{j}.polyH = dots{j}.polyH * ...
           ( H(1:bsize,1:bsize) - Shifts(1)*eye(bsize) ) + ...
           polyH * ...
           ( (dots{j}.H)(1:bsize,1:bsize) - (dots{j}.seed)(1)*eye(bsize) );
          end
          polyH *= H(1:bsize,1:bsize) - Shifts(1)*eye(bsize);
          for j=1:s
            (dots{j}.seed)(1) = [];
          end
          Shifts(1) = [];
        endfor

        for j=1:s
          dots{j}.v = (dots{j}.polyH)(:,1);
        end

        #create house vector
        [H, dots] = hhStepR(H,1,bsize,dots,sym,polyH(:,1));
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
        [H, dots] = hhStepB(H,tempStIdx,tempEndIdx,dots,sym);
        tempEndIdx = tempStIdx;
      endfor
      
      #add shift to bottom if still need to add shifts
      #The +1 is because the bulge will contain bsize - 1 shifts
      bsize = min( _MAXBULGESIZE,length(Shifts))+1;
      if(!isempty(Shifts) && (isempty(bstIdx) || bstIdx(end) + 1 > bsize ) )
        #calculate necessary polynomial bits
        for j=1:s
          dots{j}.polyH = zeros(bsize);
        end
        polyH = eye(bsize);
        #explicitly form the polynomial of shifts to shove in.
        #Note that it should be 1:bsize-1 to account for fill in
        for i=1:bsize-1
          for j=1:s
            dots{j}.polyH = dots{j}.polyH*...
              (H(n-bsize+1:n,n-bsize+1:n) - Shifts(end)*eye(bsize)) + ...
              polyH*...
              ((dots{j}.H)(n-bsize+1:n,n-bsize+1:n) - (dots{j}.seed)(end)*eye(bsize));
          end
          polyH *= H(n-bsize+1:n,n-bsize+1:n) - Shifts(end)*eye(bsize);
          for j=1:s
            (dots{j}.seed)(end) = []; 
          end
          Shifts(end) = [];
        endfor
        
        for j=1:s
          dots{j}.v = (dots{j}.polyH)(end,:);
        end
        
        [H, dots] = hhStepB(H,0,bsize-1,dots,sym,polyH(end,:));
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
        print(sprintf('../impStepPlts/impstep%03d.png',++pltNum),'-dpng');
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
  
  if(toSpike)
    INFO.HPreSchur = H;
    [H, dots, INFO] = schurComp(H, dots, spSt, spEnd, INFO, _RELTOL);
  
    if(toplt)
      pltMat(H);
      if(toprt)
        print(sprintf('../impStepPlts/impstep%03d.png',++pltNum),'-dpng');
      else
        pause(_PAUSELEN);
      endif
    endif

    #Clean up back to hessenberg form
    if(middle == 'm')
      #if enough zeros in spikes tips, deflate
      maxRow = max(find(H(:,spSt-1)));
      #check if its a complex conjugate pair
      if(abs(H(maxRow+1,maxRow)) > sqrt(eps))
        maxRow += 1;%make sure R encompasses this
      end
      minCol = min(find(H(spEnd+1,:)));
      #if it is splittable
      if(maxRow < minCol)
        #split
        Htop = H(1:maxRow,1:maxRow);
        Hbot = H(maxRow+1:end, maxRow+1:end);
        for i=1:s
          dots{i}.Htop = dots{i}.H(1:maxRow,1:maxRow);
          dots{i}.Hbot = dots{i}.H(maxRow+1:end, maxRow+1:end);
        end
        [Htop, dots, egTop] = cleanBot(Htop,dots,spSt-1,sym,toplt,toprt,INFO,_EXTRASPIKE,_RELTOL);
        [Hbot, dots, egBot] = cleanTop(Hbot,dots,spEnd+1 - maxRow,sym,toplt,toprt,INFO,_EXTRASPIKE,_RELTOL);
        #TODO: reorganize this information for use by end user
      else
        #if not enough zeros in spikes, record them and their derivates
        for i=1:s
          shiftDerivs(i,:) = [ (dots{i}.H)(spSt:spEnd,spSt-1)' (dots{i}.H)(spEnd+1,spSt:spEnd) ]; 
        end
        spikes = [ H(spSt:spEnd,spSt-1)' H(spEnd+1,spSt:spEnd) ]; 
      end
    elseif(middle == 't')#if pushed to top 
        #TODO: reorganize this information for use by end user
      [Htop, dots, egTop] = cleanTop(H,dots,spEnd+1,sym,toplt,toprt,INFO,_EXTRASPIKE,_RELTOL);
      for i=1:s
        shiftDerivs(i,:) = [ (dots{i}.H)(spEnd+1,spSt:spEnd) ];
      end
      spikes = [ H(spEnd+1,spSt:spEnd) ];
    elseif(middle == 'b') #if pushed to bottom
        #TODO: reorganize this information for use by end user
      [H, dots, spIdx, egs] = cleanBot(H,dots,spSt-1,sym,toplt,toprt,INFO,_EXTRASPIKE,_RELTOL);
      for i=1:s
        shiftDerivs(i,:) = [ (dots{i}.H)(spIdx + 1:end,spIdx)' ];
      end
      spikes = [ H(spIdx+1:end,spIdx)' ];
      [H,dots] = finishBot(H,dots,spIdx,sym,toplt,toprt);
    end

  else#not creating spikes, drive bulge off edge
    if(middle == 't')
      [H,dots] = finishTop(H,dots,n - bstIdx(end) + 1,sym,toplt,toprt);
    elseif(middle == 'b')
      [H,dots] = finishBot(H,dots,tstIdx(end),sym,toplt,toprt);
    end
  end
  INFO.dots = dots;

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
end#function

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
end#function

#perform a householder transformation, pushing a bulge rightwards
function [H, dots] = hhStepR(H,stIdx, eIdx, dots,sym,v=[])

  s = length(dots);

  #create house vector
  if(length(v) ~= eIdx - stIdx + 1)#if a vector is not provided
    for i=1:s
      dots{i}.v = (dots{i}.H)(stIdx:eIdx, stIdx-1);
    end
    v = H(stIdx:eIdx, stIdx-1);
  end
  for i=1:s
    (dots{i}.v)(1) += sgn(v(1))*v'*dots{i}.v / sqrt(v'*v);
  end
  v(1) += sgn(v(1))*sqrt(v'*v);
  for i=1:s
    dots{i}.v -= (v'*dots{i}.v) * v / (v'*v);
    dots{i}.v /= sqrt(v'*v);
  end
  v /= sqrt(v'*v);#normalize house vector

  #apply householder transformation to the right bits
  for i=1:s
    #These must come in this order (or at least the first one)
    dots{i}.H(stIdx:eIdx,:) -= v*((2*v')*dots{i}.H(stIdx:eIdx,:)) + ...
                               v*((2*dots{i}.v')*H(stIdx:eIdx,:)) + ...
                               dots{i}.v*((2*v')*H(stIdx:eIdx,:));
  end
  H(stIdx:eIdx,:) -= v*((2*v')*H(stIdx:eIdx,:));
  maxEff = min(length(H),eIdx + 1);
  for i=1:s
    #These must come in this order
    dots{i}.H(1:maxEff,stIdx:eIdx) -= (dots{i}.H(1:maxEff,stIdx:eIdx)*(2*v))*v' + ...
                                      (H(1:maxEff,stIdx:eIdx)*(2*dots{i}.v))*v' + ...
                                      (H(1:maxEff,stIdx:eIdx)*(2*v))*dots{i}.v';
  end
  H(1:maxEff,stIdx:eIdx) -= (H(1:maxEff,stIdx:eIdx)*(2*v))*v';
  if(stIdx > 1)
    H(stIdx+1:maxEff,stIdx-1) = 0;#zeros everything out for exactness
    if(sym)#if symmetric
      H(stIdx-1,stIdx+1:maxEff) = 0;#zeros everything out for exactness
    end

    for i=1:s
      (dots{i}.H)(stIdx+1:maxEff,stIdx-1) = 0;#zeros everything out for exactness
      if(sym)#if symmetric
        (dots{i}.H)(stIdx-1,stIdx+1:maxEff) = 0;#zeros everything out for exactness
      end
    end
  end
end

#perform a householder transformation, pushing a bulge from bottom
function [H, dots] = hhStepB(H,stIdx, eIdx,dots,sym,v=[])
  n = length(H);
  s = length(dots);
 
  #create house vector
  if(length(v) ~= eIdx - stIdx + 1)#if a vector is not provided
    for i=1:s
      dots{i}.v = (dots{i}.H)(n - stIdx + 1, n - eIdx:n - stIdx);
    end
    v = H(n - stIdx + 1, n - eIdx:n - stIdx);
  end
  for i=1:s
    (dots{i}.v)(end) += sgn(v(end))*dots{i}.v*v'/sqrt(v*v');
  end
  v(end) += sgn(v(end))*sqrt(v*v');
  for i=1:s
    dots{i}.v -= (dots{i}.v*v') * v/(v*v');
    dots{i}.v /= sqrt(v*v');
  end
  v /= sqrt(v*v');#normalize house vector
  
  #apply householder transformation to the right bits
  minEff = max(1,n - eIdx - 1);
  for i=1:s
    dots{i}.H(n - eIdx:n - stIdx, minEff:n) -= v'*((2*v)*dots{i}.H(n - eIdx:n - stIdx,minEff:n)) + ...
                                               v'*((2*dots{i}.v)*H(n - eIdx:n - stIdx,minEff:n)) + ...
                                               dots{i}.v'*((2*v)*H(n - eIdx:n - stIdx,minEff:n));
  end
  H(n - eIdx:n - stIdx, minEff:n) -= v'*((2*v)*H(n - eIdx:n - stIdx,minEff:n));
  for i=1:s
    #The order here matters
    dots{i}.H(:,n - eIdx:n - stIdx) -= (dots{i}.H(:,n - eIdx:n - stIdx)*(2*v'))*v + ...
                                       (H(:,n - eIdx:n - stIdx)*(2*dots{i}.v'))*v + ...
                                       (H(:,n - eIdx:n - stIdx)*(2*v'))*dots{i}.v;
  end
  H(:,n - eIdx:n - stIdx) -= (H(:,n - eIdx:n - stIdx)*(2*v'))*v;
  if(stIdx > 0)
    H(n - stIdx + 1,minEff: n - stIdx-1) = 0;#zeros everything out for exactness
    if(sym)#if symmetric
      H(minEff: n - stIdx-1,n - stIdx + 1) = 0;#zeros everything out for exactness
    end

    for i=1:s
      (dots{i}.H)(n - stIdx + 1,minEff: n - stIdx-1) = 0;#zeros everything out for exactness
      if(sym)#if symmetric
        (dots{i}.H)(minEff: n - stIdx-1,n - stIdx + 1) = 0;#zeros everything out for exactness
      end
    end
  end
end

#perform the schur decomposition
function [H,dots,INFO] = schurComp(H,dots,spSt,spEnd,INFO,_RELTOL)
  
  n = length(H);
  s = length(dots);
  INFO.nz = 0;

  #actually perform the schur decomposion (AU = UT)
  [spRot, H(spSt:spEnd,spSt:spEnd)] = schur(H(spSt:spEnd,spSt:spEnd));

  middle = 'm';
  if(spSt == 1)
    middle = 't';
  elseif(spEnd == length(H))
    middle = 'b';
  end
  #sort the decomposition to achieve good spikes
  [Q, T, ap] = swapSchur(spRot, H(spSt:spEnd,spSt:spEnd), middle, _RELTOL);

  if(max(ap) > 1)
    warning('schur reordering may have performed badly,so did not reorder')
  else
    spRot = Q;
    H(spSt:spEnd,spSt:spEnd) = T;
  end

  for i=1:s
    [dots{i}.spRot, (dots{i}.H)(spSt:spEnd,spSt:spEnd)] = ...
      schurAd((dots{i}.H)(spSt:spEnd,spSt:spEnd),spRot, H(spSt:spEnd,spSt:spEnd));
  end
  if(spSt > 1)%apply it to the left/top portion
    for i=1:s
      (dots{i}.H)(1:spSt-1,spSt:spEnd) = ...
        (dots{i}.H)(1:spSt-1,spSt:spEnd) * spRot + H(1:spSt-1,spSt:spEnd) * dots{i}.spRot;
    end
    H(1:spSt-1,spSt:spEnd) = H(1:spSt-1,spSt:spEnd) * spRot;
  end
  if(spEnd < n)%apply it to the right/bottom portion
    for i=1:s
      (dots{i}.H)(spSt:spEnd,spEnd+1:n) = ...
        dots{i}.spRot'*H(spSt:spEnd,spEnd+1:n) + spRot'*(dots{i}.H)(spSt:spEnd,spEnd+1:n);
    end
    H(spSt:spEnd,spEnd+1:n) = spRot'*H(spSt:spEnd,spEnd+1:n);
  end

  #form spikes
  vSpike = spRot'(:,1);
  hSpike = spRot(end,:);
  for i=1:s
    dots{i}.vSpike = (dots{i}.spRot')(:,1);
    dots{i}.hSpike = (dots{i}.spRot)(end,:);
  end
  #count hard zeros
  vsz = abs(vSpike) < _RELTOL;
  hsz = abs(hSpike) < _RELTOL;
  #enforce hard zeros
  #This is non-smooth, so not performed when differentiating
  if(s == 0)
    vSpike(vsz) = 0;
    hSpike(hsz) = 0;
  end
  
  if(spEnd < n)%apply it to the right/bottom portion
    for i=1:s
      (dots{i}.H)(spEnd+1,spSt:spEnd) = ...
        (dots{i}.H)(spEnd+1,spEnd) * hSpike + H(spEnd+1,spEnd) * dots{i}.hSpike;;
    end
    H(spEnd+1,spSt:spEnd) = H(spEnd+1,spEnd) * hSpike;
    INFO.nz += sum(hsz);
    INFO.hspIdx = spEnd + 1;
  end
  if(spSt > 1)%apply it to the left/top portion
    for i=1:s
      (dots{i}.H)(spSt:spEnd,spSt-1) = ...
        (dots{i}.H)(spSt,spSt-1) * vSpike + H(spSt,spSt-1) * dots{i}.vSpike;
    end
    H(spSt:spEnd,spSt-1) = H(spSt,spSt-1) * vSpike;
    INFO.nz += sum(vsz);
    INFO.vspIdx = spSt - 1;
  end

end

#finishes sweep to the bottom of a matrix
function [H, dots, spIdx, eg] = cleanBot(H, dots, spIdx, sym, toplt, toprt,INFO, _EXTRASPIKE, _RELTOL)
  s = length(dots);
  n = length(H);
  defIdx = max(find(H(:,spIdx)))+1;#topmost zero
  spLength = n - spIdx;
  if(defIdx <= n && abs(H(defIdx,defIdx-1)) > _RELTOL)
    #this represents the bottom portion of a complex
    #conjugate 2x2 block, and the top half of it is not
    #deflatable, so skip it
    defIdx++;
  end
  B = [];
  tryAgain = false;
  if(defIdx <= n)
    B = H(defIdx:n,defIdx:n);
    tryAgain = true;
  end
  H = H(1:defIdx-1, 1:defIdx-1);
  for i=1:s
    dots{i}.H = dots{i}.H(1:defIdx-1, 1:defIdx-1);
  end
  eg = [];
  if( ~isempty(B) )#if deflatable
    #grab the eigenvalues from the deflated matrix
    while(~isempty(B))
      if(length(B) == 1)
        eg = [eg,B];
        B = [];
      else
        if(abs(B(2,1)) > _RELTOL)
          #extract complex cojugate pair
          rp = B(1,1) + B(2,2);
          rp *= .5;
          ip = rp^2 - B(1,1)*B(2,2) + B(1,2)*B(2,1);
          ip = sqrt(ip);
          eg = [eg, rp + ip, rp - ip];
          B = B(3:end,3:end);
        else
          #extract the eigenvalue
          eg = [eg,B(1,1)];
          B = B(2:end,2:end);
        end
      end
    end
  end


  if(tryAgain)
    [H, dots] = finishBot(H,dots,spIdx+1,sym,toplt,toprt);
    #respike
    n = length(H);
    spSt = n - spLength;#spike will form 1 left of here
    spEnd = n;
    [H,dots,INFO] = schurComp(H,dots,spSt,spEnd,INFO,_RELTOL);

    spIdx = spSt - 1;
    #perform this again
    [H, dots, spIdx, eg2] = cleanBot(H,dots,spIdx,sym,toplt,toprt,INFO,_EXTRASPIKE,_RELTOL);
    eg = [eg, eg2];
  end

end

#finishes a sweep to the top of a matrix
function [H, dots, eg] = cleanTop(H, dots, spIdx, sym, toplt, toprt, INFO,_EXTRASPIKE, _RELTOL)
  s = length(dots);
  defIdx = min(find(H(spIdx,:))) - 1;#rightmost zero
  spLength = spIdx;
  if(defIdx >= 1 && abs(H(defIdx+1,defIdx)) > _RELTOL)
    #this represents the bottom portion of a complex
    #conjugate 2x2 block, and the top half of it is not
    #deflatable, so skip it
    defIdx--;
  end
  T = [];
  tryAgain = false;
  if(defIdx >= 1)
    T = H(1:defIdx,1:defIdx);
    tryAgain = true;
  end
  H = H(defIdx+1:end, defIdx+1:end);
  for i=1:s
    dots{i}.H = dots{i}.H(defIdx+1:end, defIdx+1:end);
  end
  eg = [];
  if( ~isempty(T) )#if deflatable
    #grab the eigenvalues from the deflated matrix
    while(~isempty(T))
      if(length(T) == 1)
        eg = [eg,T];
        T = [];
      else
        if(abs(T(2,1)) > _RELTOL)
          #extract complex cojugate pair
          rp = T(1,1) + T(2,2);
          rp *= .5;
          ip = rp^2 - T(1,1)*T(2,2) + T(1,2)*T(2,1);
          ip = sqrt(ip);
          eg = [eg, rp + ip, rp - ip];
          T = T(3:end,3:end);
        else
          #extract the eigenvalue
          eg = [eg,T(1,1)];
          T = T(2:end,2:end);
        end
      end
    end
  end

  [H, dots] = finishTop(H,dots,spIdx,sym,toplt,toprt);

  if(tryAgain)
    #respike
    spSt = 1;
    spEnd = spLength;
    [H, dots, INFO] = schurComp(H, dots, spSt, spEnd, INFO, _RELTOL);

    #perform this again
    [H,dots,eg2] = cleanTop(H,dots,spIdx,sym,toplt,toprt);
    eg = [eg, eg2];
  end

end

function [H, dots] = finishTop(H,dots,stIdx, sym, toplt, toprt)
  #reduce H to hessenberg form
  n = length(H);
  pltNum = 0;
  for toHess = stIdx:-1:3
    #The indices look funny because they were desinged to take
    #input from the trains and we don't need trains here
    [H, dots] = hhStepB(H,n+1-toHess,n-1,dots,sym);
    if(toplt)
      pltMat(H);
      if(toprt)
        print(sprintf('../impStepPlts/impstepFinish%03d.png',++pltNum),'-dpng');
      else
        pause(0.0);
      endif
    endif
  end
end

function [H, dots] = finishBot(H,dots,stIdx,sym,toplt,toprt)
  #reduce H to hessenberg form
  n = length(H);
  pltNum = 0;
  for toHess = stIdx:length(H)
    [H, dots] = hhStepR(H,toHess,n,dots,sym);
    if(toplt)
      pltMat(H);
      if(toprt)
        print(sprintf('../impStepPlts/impstepFinish%03d.png',++pltNum),'-dpng');
      else
        pause(0.0);
      endif
    endif
  end
end
