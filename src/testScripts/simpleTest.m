#Hermitian test
eVals = 1 - 2*rand(100,1);
[Q,~] = qr(rand(100));
A = hess(Q' * diag(eVals) * Q);

printf('testing bottom deflation with 15 eigenvalue shifts\n');
[H,INFO] = trainBust(A,eVals(1:15),false,false,'b'); 
pltMat(H); 
printf('The spike had %d zeros\n',INFO.nz);
printf('max difference between eigenvalues pre and post train run  = %g\n',max(abs(sort(eig(H)) - sort(eVals))));
printf('hit a key to continue\n\n');
pause();

printf('testing top deflation with 15 eigenvalue shifts\n');
[H,INFO] = trainBust(A,eVals(1:15),false,false,'t'); 
pltMat(H); 
printf('The spike had %d zeros\n',INFO.nz);
printf('max difference between eigenvalues pre and post train run  = %g\n',max(abs(sort(eig(H)) - sort(eVals))));
printf('hit a key to continue\n\n');
pause();

printf('testing middle deflation with 15 eigenvalue shifts\n');
[H,INFO] = trainBust(A,eVals(1:15),false,false,'m'); 
pltMat(H); 
printf('The spike had %d zeros\n',INFO.nz);
printf('max difference between eigenvalues pre and post train run  = %g\n',max(abs(sort(eig(H)) - sort(eVals))));
printf('hit a key to continue\n\n');
pause();

#Non-hermitian test
eVals = 1 - 2*rand(100,1);
Q = rand(100);
A = hess(inv(Q) * diag(eVals) * Q);

printf('testing bottom deflation with 15 eigenvalue shifts\n');
[H,INFO] = trainBust(A,eVals(1:15),false,false,'b'); 
pltMat(H); 
printf('The spike had %d zeros\n',INFO.nz);
printf('max difference between eigenvalues pre and post train run  = %g\n',max(abs(sort(eig(H)) - sort(eVals))));
printf('hit a key to continue\n\n');
pause();

printf('testing top deflation with 15 eigenvalue shifts\n');
[H,INFO] = trainBust(A,eVals(1:15),false,false,'t'); 
pltMat(H); 
printf('The spike had %d zeros\n',INFO.nz);
printf('max difference between eigenvalues pre and post train run  = %g\n',max(abs(sort(eig(H)) - sort(eVals))));
printf('hit a key to continue\n\n');
pause();

printf('testing middle deflation with 15 eigenvalue shifts\n');
[H,INFO] = trainBust(A,eVals(1:15),false,false,'m'); 
pltMat(H); 
printf('The spike had %d zeros\n',INFO.nz);
printf('max difference between eigenvalues pre and post train run  = %g\n',max(abs(sort(eig(H)) - sort(eVals))));
printf('hit a key to continue\n\n');
pause();
