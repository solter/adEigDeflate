n = 1000;
vals = nan(n,2);
for tomid = [true, false]
  dirt = '';
  if(tomid)
    dirt = '../zeroMidPlots';
  else
    dirt = '../zeroBotPlots';
  endif
  mkdir(dirt);
  er = nan;
  do#do while still getting holes in spike
    for i=1:n#for each sample
      A = 1 - 2*rand(140);
      ev = 1 - 2*rand(140,1);
      A = inv(A)*diag(ev)*A;
      A = hess(A);
     
      sft = zeros(14,1);
      #shifts 
      if(!isnan(er))#varied from evals
        sft = 1 - 2*rand(14,1);
        sft *= 10^er;
        sft += 1;
        sft = sft .* ev(1:14);
      else#which are evals
        sft = ev(1:14);
      endif


      [H,info] = trainBust(A, sft,false,false,tomid);
      vals(i,1) = info.nz;
      #vals(i,2) = norm(sort(real(eig(H))) - sort(ev));
      if(mod(i,floor(n/10)) == 0)
        printf("\n%.1f%% finished wih er = %.1f\n",100*i/n,er);
        fflush(stdout);
      elseif(mod(i,10) == 0)
        printf(".");
        fflush(stdout);
      endif
    endfor
    figure();
    hist(vals(:,1),min(vals(:,1)):max(vals(:,1))+1,100);
    ylabel("%")
    if(!isnan(er))
      title(sprintf("Num 0's with %.1g variation in evals", 10^(er)))
      print(sprintf('%s/1e%.1f.png',dirt,er))
      close all
    else
      title("Num 0's with exact eigenvalues")  
      print(sprintf('%s/exact.png',dirt))
      close all
      er = -9;
    endif
    printf("\nDone with error = %d\n",er);
    fflush(stdout);
    er++;
  until(sum(vals(:,1)) == 0)
  printf("\n\nDone with 1 cycle of tomid\n");
  fflush(stdout)
endfor
