n = 1000;
mkdir('../zeroPlots')
vals = nan(n,2);
for er = [nan, -4:.5:-1]
  parfor i=1:n
    A = 1 - 2*rand(140);
    ev = 1 - 2*rand(140,1);
    A = inv(A)*diag(ev)*A;
    A = hess(A);
    
    #shifts 
    if(er != nan)#varied from evals
      sft = 1 - 2*rand(14,1);
      sft *= 10^er;
      sft += 1;
      sft = sft .* ev(1:14);
    else#which are evals
      sft = ev(1:14);
    endif


    [H,info] = trainBust(A, sft,false,false,false);
    vals(i,1) = info.nz;
    #vals(i,2) = norm(sort(real(eig(H))) - sort(ev));
    if(mod(i,100) == 0)
      printf("\n%d%% finished wih er = %.1f\n",i/10,er);
      fflush(stdout);
    elseif(mod(i,10) == 0)
      printf(".");
      fflush(stdout);
    endif
  endparfor
  figure();
  hist(vals(:,1),min(vals(:,1)):max(vals(:,1)),100);
  ylabel("%")
  if(er != nan)
    title(sprintf("Num 0's with %.2f %% in evals", 100*10^(-er)))
    print('../zeroPlots/1e%.1f.png',er)
  else
    title("Num 0's with exact eigenvalues")  
    print('../zeroPlots/exact.png')
  endif
endfor
