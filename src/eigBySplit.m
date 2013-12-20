##-*- texinfo -*-
##@deftypefn{Function File} {@var{EVal}=} eigBySplit(@var{A})
##@deftypefnx{Function File} {@var{EVal}=} eigBySplit(@var{A}, @var{isHess})
##@deftypefnx{Function File} {@var{EVal}=} eigBySplit(@var{A}, @var{isHess}, @var{shifts})
##
##Solves for the eigenvalues of A.
##
##This is done through a reduction of A to hessenberg form.
##Then shifts are chosen through an optimization strategy
##which allows for the deflation of the matrix.
##This is done recursively until the matrix is small.
##
##Inputs:@*
##  @var{A} - matrix to find evals for @*
##  @var{isHess} - if true, will not bother converting it to hessenberg form internally
##  @var{shifts} - optional choice for starting shifts.
##    Will be split in half with the first half applied to the top
##    Should be less than ~12 to avoid shift blurring.@*
##
##Outputs:@*
##  @var{EVal} - A list of the eigenvalues.@* 
## @end deftypefn
function [EVal] = eigBySplit(A, isHess = false, shifts = [])
  if(!isHess)
    H = hess(A);
  else
    H = A;
  endif

  #TODO: make relative tolerance
  tol = 1e-2;

  s = 5;
  if(isempty(shifts))
    shifts = real(eig(H(end-2*s:end,end-2*s:end)));
  else
    s = floor(s/2);
  endif

  #try one step
  #create spikes
  [H, sp1, sp2] = colAndBust(H,shifts(1:s),shifts(s+1:end));

  #extract spike values
  spike = [H(sp1+2:sp2-1,sp1); H(sp2, sp1+1:sp2-2)']
    
  #+++++optimize shifts+++++
  while(0 && sum(abs(spike) < tol) < length(shifts)) #indicates good shifts
    #TODO: update shift estimate

    #create spikes
    [H, sp1, sp2] = colAndBust(H,shifts(1:s),shifts(s+1:end));

    #extract spike values
    spike = [H(sp1+2:sp2-1,sp1); H(sp2, sp1+1:sp2-2)'];
    
  endwhile
  EVal.sp = H;  
  #recollapse to hessenberg
  H = hess(H);

  EVal.h = H;
endfunction
