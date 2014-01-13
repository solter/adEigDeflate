##-*- texinfo -*-
##@deftypefn{Function File} {[@var{H}, @var{spSt}, @var{spEnd}, @var{pltNum}] =} agErlyDef(@var{A},@var{topshifts})
##@deftypefnx{Function File} {[@var{H}, @var{spSt}, @var{spEnd}, @var{pltNum}] =} agErlyDef(@var{A},@var{topshifts},@var{toplt})
##@deftypefnx{Function File} {[@var{H}, @var{spSt}, @var{spEnd}, @var{pltNum}] =} agErlyDef(@var{A},@var{topshifts},@var{toplt},@var{toprt})
##@deftypefnx{Function File} {[@var{H}, @var{spSt}, @var{spEnd}, @var{pltNum}] =} agErlyDef(@var{A},@var{topshifts},@var{toplt},@var{toprt},@var{spike})
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
##  @var{spike} - optional argument. If it is false, then the bulges
##    will be chased off the bottom rather than creating a spike
##
##Outputs:@*
##  @var{H} - The matrix A after running 2 bulges into it's center and
##    creating the spikes@*
## @end deftypefn
function [H, pltNum] = agErlyDef(H, Shifts, toplt = false, toprt = false, spike = true)

#TODO: actually deflate matrix, so far only finds deflation points while pushing bulges
#need to do deflation logic for spike

#TODO 2: implement derivatives... use in deflation logic

  _PAUSELEN = 0.1;#pause time when playing real time
  _MAXBULGESIZE = 4;#maximum bulge size
  _EXTRASPIKE = .5*length(Shifts);#extra size of spike beyond shifts
  _RELTOL = 1e-6;
  _ABSTOL = 1e-14;

  toprt = toprt && tplt;%must be plotting to print

  if(toprt)
    mkdir('../impStepPlts');
    pltNum = 0;
  else
    toprt = false;
  end#if

  stIdx = [];#bulge starts
  defIx = [];#list of points to deflate at
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

      if(++bsize >= _MAXBULGESIZE || isempty(Shifts))
        stIdx = [stIdx 1];
        bsize = 0;
      endif
    elseif( stIdx(end) > 3 )
      #check for deflatable point
      potIdx = stIdx(end) - 2;
      #if entry is sufficientyly small
      if ( abs( H(potIdx+1,potIdx)) < _RELTOL*( abs( H(potIdx,potIdx) )  + abs( H(potIdx + 1, potIdx + 1) ) ) || abs( H(potIdx+1,potIdx) ) < _ABSTOL)
        defIdx = [potIdx defIdx];
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
  until(endIdx + 1 >= length(H))#if it touches the bottom

  if(spike)

    #create spike
    spSt = stIdx(end) - _EXTRASPIKE;#index of block
    [spRot,~] = schur(H(spSt:end,spSt:end));
    H(:,spSt:end) = H(:,spSt:end)*spRot;
    H(spSt:end,(spSt-1):end) = spRot'*H(spSt:end,(spSt-1):end);

    if(toplt)
      pltMat(H);
      if(toprt)
        print(sprintf('../impStepPlts/impstep%03d.png',++pltNum));
      else
        pause(4*_PAUSELEN+.1);
      endif
    end#if

    #TODO: Deflate

    #chase spike off matrix
    stIdx = spSt;
    do

      v = H(stIdx:end, stIdx-1);
      v(1) += sgn(v(1))*sqrt(v'*v);
      v /= sqrt(v'*v);#normalize house vector
      #apply householder transformation to the right bits
      H(stIdx:end,:) -= v*((2*v')*H(stIdx:end,:));
      H(1:end,stIdx:end) -= (H(1:end,stIdx:end)*(2*v))*v';
      H(stIdx+1:end,stIdx-1) = 0;#zeros everything out for exactness

      stIdx++;

      if(toplt)
        pltMat(H);
        if(toprt)
          print(sprintf('../impStepPlts/impstep%03d.png',++pltNum));
        else
          pause(_PAUSELEN);
        endif
      end#if

    until(stIdx >= length(H))

  else#chase bulge off end

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
          pause(_PAUSELEN);
        endif
      end#if

    until(isempty(stIdx))

  endif%done with spiking branching

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
