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
##    specified by @var{shiftSeed}. 
##  @var{INFO} - a structure with INFOrmation. Fields:
##    vspIdx -> vertical spike's column index
##    hspIdx -> horizontal spike's row index (only for middle)
##    nz -> count of zeros in spikes
## @end deftypefn
function [H, shiftDerivs, INFO] = trainBust(H, Shifts, shiftSeed, toplt = false, toprt = false, middle = 'm')

  #NOTE: all the derivatives MUST be updated before the quantity is

  #TODO: fix the (I-vv')H(I-vv') derivatives. There's an issue in that
  #we really only form the small non-zero v, which fails to accurately
  #recreate the last 2 terms as currently coded

  #TODO: actually deflate matrix, so far only finds deflation points while pushing bulges
  #need to do deflation logic for spike

  _PAUSELEN = 0.0;#pause time when playing real time
  _MAXBULGESIZE = 5;#maximum bulge size
  _EXTRASPIKE = 0;#floor(.5*length(Shifts));#extra size of spike beyond shifts
  _RELTOL = 1e-6;
  n = length(H);#convenience shorthand
  toprt = toprt && toplt;#must be plotting to print
  INFO = struct();
  sym = norm(H - H','inf') < sqrt(eps);
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
        polyH = eye(bsize);
        for j=1:s
          dots{j}.polyH = zeros(bsize);
        end
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
        polyH = eye(bsize);
        for j=1:s
          dots{j}.polyH = zeros(bsize);
        end
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
    
  [H, dots, INFO] = schurComp(H, dots, spSt, spEnd, INFO, _RELTOL);
  
  #actually extract the appropriate spikes and their shifts
  shiftDerivs = [];
  for i=1:s
    if(middle == 'm')
      shiftDerivs(i,:) = [ (dots{i}.H)(spSt:spEnd,spSt-1)' (dots{i}.H)(spEnd+1,spSt:spEnd) ]; 
    elseif(middle == 'b')
      shiftDerivs(i,:) = [ (dots{i}.H)(spSt:spEnd,spSt-1)' ];
    elseif(middle == 't')
      shiftDerivs(i,:) = [ (dots{i}.H)(spEnd+1,spSt:spEnd) ];
    end
  end


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
function [H, dots] = hhStepR(H,stIdx, eIdx, dots,sym,v=[])

  s = length(dots);
  one = 1:stIdx-1;
  two = stIdx:eIdx;
  three = eIdx+1:length(H);

  #create house vector
  if(length(v) ~= eIdx - stIdx + 1)#if a vector is not provided
    for i=1:s
      dots{i}.v = (dots{i}.H)(two, stIdx-1);
    end
    v = H(two, stIdx-1);
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
  maxEff = min(length(H),eIdx + 1);
  for i=1:s
    #This is structures as it is because
    #the 2 steps left and right multiplying by 2vv' is really
    #performing (I-2vv')H(I-2vv'). Thus the derivative has lots of terms.
    #
    #The following algorithm was generated by considering matrix
    #as a 3x3 block matrix, where the central block is indices (stIdx:eIdx,stIdx:eIdx) = (two,two).
    #Then vv'H only affects the second row, and Hvv' affects only the second column.
    #Expanding these led to the following algorithm (noting that all vv' are symmetric)
    #
    #Note that the corners are not effected

    #This could be much more efficient if the hessenberg + bulge structure is exploited
    #All the extra parenthesis ensure no matrix matrix products are ever performed, at worst it's
    #a matrix of the form a*b' formed

    #dots.H32  = dots.H32*(2vv'-I) + 2*H32*(dots.v*v' + v*dots.v') - 4*H32_updated
    #Note that the updated needs to happen after H is updated
    (dots{i}.H)(three,two) = (dots{i}.H)(three,two)*2*v*v' - (dots{i}.H)(three,two) + 
                            2*(H(three,two) *dots{i}.v)*v' + 2*(H(three,two)*v)*dots{i}.v');
    if(sym)
      (dots{i}.H)(two,three) = (dots{i}.H)(three,two)';
    else
      (dots{i}.H)(one,two) = (dots{i}.H)(one,two)*2*v*v' - (dots{i}.H)(one,two) + 
                              2*(H(one,two) *dots{i}.v)*v' + 2*(H(one,two)*v)*dots{i}.v');
    end
    #dots.H21 = (2vv'-I)*dots.H21 + 2*(dots.v*v' + v*dots.v')*H21 - 4*H21_updated
    #Note that the updated needs to happen after H is updated
    (dots{i}.H)(two,one) = 2*v*v'*(dots{i}.H)(two,one) - (dots{i}.H)(two,one) + 
                         2*dots{i}.v*(v'*H(two,one)) + 2*v*(dots{i}.v'*H(two,one));
    if(sym)
      (dots{i}.H)(one,two) = (dots{i}.H)(two,one)'; 
    else
      (dots{i}.H)(two,three) = 2*v*v'*(dots{i}.H)(two,three) - (dots{i}.H)(two,three) + 
                           2*dots{i}.v*(v'*H(two,three)) + 2*v*(dots{i}.v'*H(two,three));
    end
    #dots.H22 = (I-2vv')*dots.H22*(I-2vv') - 2(I - 2vv')H22*(dots.v*v' + v*dots.v') - 2(dots.v*v' + v*dots.v')*H22*(I-2vv')
    #these are expanded out, and the last two terms are combined if its symmetric
    (dots{i}.H)(two,two) = (dots{i}.H)(two,two) - 2*(v*(v'*(dots{i}.H)(two,two)) + ((dots{i}.H)(two,two)*v)*v') + 4*v*(v'*((dots{i}.H)(two,two)*v))*v';
    if(sym)
      (dots{i}.H)(two,two) -= 4*( dots{i}.v*(v'*H(two,tow)) + v*(dots{i}.v*H(two,two)) + 2*dots{i}.v*(v'*H(two,two)*v)*v' + 2*v*(dots{v}.v'*H(two,two)*v)*v' );
    else
      (dots{i}.H)(two,two) -= 2*( dots{i}.v*(v'*H(two,tow)) + v*(dots{i}.v*H(two,two)) + 2*dots{i}.v*(v'*H(two,two)*v)*v' + 2*v*(dots{v}.v'*H(two,two)*v)*v' );
      #Since everything but H is symmetric, the following is equivalent
      (dots{i}.H)(two,two) -= 2*( dots{i}.v*(v'*H(two,tow)') + v*(dots{i}.v*H(two,two)') + 2*dots{i}.v*(v'*H(two,two)'*v)*v' + 2*v*(dots{v}.v'*H(two,two)'*v)*v' )';
    end
  end
  H(two,:) -= v*((2*v')*H(two,:));
  H(1:maxEff,two) -= (H(1:maxEff,two)*(2*v))*v';
  for i=1:s
    #The updated parts
    #TODO
  end
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
    #This looks kinda funny for the following reasons
    #1) the 2 steps left and right multiplying by 2vv' is really
    #performing (I-2vv')H(I-2vv').
    #When you expand this out you get the following:

    #The next 3 can come in any order
    (dots{i}.H)(n - eIdx:n - stIdx, minEff:n) -= dots{i}.v'*((2*v)*H(n - eIdx:n - stIdx,minEff:n));
    (dots{i}.H)(n - eIdx:n - stIdx, minEff:n) -= v'*((2*dots{i}.v)*H(n - eIdx:n - stIdx,minEff:n));
    (dots{i}.H)(n - eIdx:n - stIdx, minEff:n) -= v'*((2*v)*(dots{i}.H)(n - eIdx:n - stIdx,minEff:n));
    
    #This one MUST come next
    (dots{i}.H)(:,n - eIdx:n - stIdx) -= ((dots{i}.H)(:,n - eIdx:n - stIdx)*(2*v'))*v;
    
    #The next 4 can come in any order
    (dots{i}.H)(:,n - eIdx:n - stIdx) -= (H(:,n - eIdx:n - stIdx)*(2*dots{i}.v'))*v;
    (dots{i}.H)(:,n - eIdx:n - stIdx) -= (H(:,n - eIdx:n - stIdx)*(2*v'))*dots{i}.v;
    (dots{i}.H)(:,n - eIdx:n - stIdx) -= v'*((4*v)*(H(:,n - eIdx:n - stIdx)*(dots{i}.v')))*v;
    (dots{i}.H)(:,n - eIdx:n - stIdx) -= v'*((4*v)*(H(:,n - eIdx:n - stIdx)*(v')))*dots{i}.v;
  end
  H(n - eIdx:n - stIdx, minEff:n) -= v'*((2*v)*H(n - eIdx:n - stIdx,minEff:n));
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

  #actually perform the schur decomposiotn
  [spRot, H(spSt:spEnd,spSt:spEnd)] = schur(H(spSt:spEnd,spSt:spEnd));
  for i=1:s
    [dots{i}.spRot, (dots{i}.H)(spSt:spEnd,spSt:spEnd)] = ...
      schurAd((dots{i}.H)(spSt:spEnd,spSt:spEnd),spRot, H(spSt:spEnd,spSt:spEnd));
  end
  if(spSt > 1)%apply it to the left/top portion
    H(1:spSt-1,spSt:spEnd) = H(1:spSt-1,spSt:spEnd) * spRot;
    for i=1:s
      (dots{i}.H)(1:spSt-1,spSt:spEnd) = ...
        (dots{i}.H)(1:spSt-1,spSt:spEnd) * spRot + H(1:spSt-1,spSt:spEnd) * dots{i}.spRot;
    end
  end
  if(spEnd < n)%apply it to the right/bottom portion
    H(spSt:spEnd,spEnd+1:n) = spRot'*H(spSt:spEnd,spEnd+1:n);
    for i=1:s
      (dots{i}.H)(spSt:spEnd,spEnd+1:n) = ...
        dots{i}.spRot'*H(spSt:spEnd,spEnd+1:n) + spRot'*(dots{i}.H)(spSt:spEnd,spEnd+1:n);
    end
  end

  #form spikes
  vSpike = spRot'(:,1);
  hSpike = spRot(end,:);
  for i=1:s
    dots{i}.vSpike = dots{i}.spRot'(:,1);
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
    H(spEnd+1,spSt:spEnd) = H(spEnd+1,spEnd) * hSpike;
    for i=1:s
      (dots{i}.H)(spEnd+1,spSt:spEnd) = ...
        (dots{i}.H)(spEnd+1,spEnd) * hSpike + H(spEnd+1,spEnd) * dots{i}.hSpike;;
    end
    INFO.nz += sum(hsz);
    INFO.hspIdx = spEnd + 1;
  end
  if(spSt > 1)%apply it to the left/top portion
    H(spSt:spEnd,spSt-1) = H(spSt,spSt-1) * vSpike;
    for i=1:s
      (dots{i}.H)(spSt:spEnd,spSt-1) = ...
        (dots{i}.H)(spSt,spSt-1) * vSpike + H(spSt,spSt-1) * dots{i}.vSpike;
    end
    INFO.nz += sum(vsz);
    INFO.vspIdx = spSt - 1;
  end

end
