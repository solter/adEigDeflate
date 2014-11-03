nzs = zeros(6,1);
niter = 200;

parfor i=1:niter
  #hermitian
  eVals = 1 - 2*rand(100,1);
  [Q,~] = qr(rand(100));
  A = hess(Q' * diag(eVals) * Q);
  [~,INFO] = trainBust(A,eVals(1:15),false,false,'b'); 
  nzs(1) += INFO.nz;
  [~,INFO] = trainBust(A,eVals(1:15),false,false,'m'); 
  nzs(2) += INFO.nz;
  [~,INFO] = trainBust(A,eVals(1:15),false,false,'t'); 
  nzs(3) += INFO.nz;

  #nonhermitian  
  eVals = 1 - 2*rand(100,1);
  Q = rand(100);
  A = hess(inv(Q) * diag(eVals) * Q);
  [~,INFO] = trainBust(A,eVals(1:15),false,false,'b'); 
  nzs(4) += INFO.nz;
  [~,INFO] = trainBust(A,eVals(1:15),false,false,'m'); 
  nzs(5) += INFO.nz;
  [~,INFO] = trainBust(A,eVals(1:15),false,false,'t'); 
  nzs(6) += INFO.nz;

  if(mod(i,niter/10) == 0)
    printf('%d tests performed\n',i);
    fflush(stdout);
  endif

endparfor

nzs /= niter;

printf('With 15/100 evals used as shifts, the number of zeros on spikes was\nbge loc \t | \t hermitian \t | \t non-herm.\n  b \t | \t %d \t | \t %d \n  m \t | \t %d \t | \t %d \n  t \t | \t %d \t | \t %d \n', nzs(1), nzs(4), nzs(2), nzs(5), nzs(3), nzs(6));

