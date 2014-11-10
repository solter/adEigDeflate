#This script is meant to test derivatives

#test the schurAd function
if(exist('tSchur','var') && tSchur == true)

  Err = [];
  ErrSort = [];
  ErrAbs = [];

  for i=1:100
    A = randn(20);
    D = rand(20);
    %orthonormalize
    [D,~] = qr(D);
    ep = 1e-8;
    [Q1,S1] = schur(A + ep*D);
    [Q2,S2] = schur(A - ep*D);
    [Q,S] = schur(A);
    [Qd,Td] = schurAd(D,Q,S);

    %because we only really care about these values:
    int1 = [Q1'(:,1)' Q1(end,:)];
    int2 = [Q2'(:,1)' Q2(end,:)];
    int = [Q'(:,1)' Q(end,:)];
    Qdint = [Qd'(:,1)' Qd(end,:)];
    [~,idxSort] = sort(int);
    [~,idxAbs] = sort(abs(int));
    
    QdCd = (int1 - int2) / (2*ep);#straight up central difference
    QdCdSort = (sort(int1) - sort(int2)) / (2*ep);#sorted central difference
    QdCdAbs = (sort(abs(int1)) - sort(abs(int2))) / (2*ep);#sorted central difference

    Err = [Err, abs(QdCd - Qdint)];
    ErrSort = [ErrSort, abs(QdCdSort - Qdint(idxSort))];
    ErrAbs = [ErrAbs,abs(QdCdAbs - Qdint(idxAbs))];
    ErrRel = [Err, abs(QdCd - Qdint)./abs(Qdint)];
    ErrSortRel = [ErrSort, abs(QdCdSort - Qdint(idxSort))./(abs(Qdint(idxSort)))];
    ErrAbsRel = [ErrAbs,abs(QdCdAbs - Qdint(idxAbs))./abs(Qdint(idxSort))];
  end

end

#tSchur2 - test analogous to the one done in ADTalkFeb4 by Dr. Struthers

if(exist('tSchur2','var') && tSchur2 == true)
  A = 2*rand(20)-1;
  Ad = 2*rand(20)-1;
  [Q,T] = schur(A);
  [Qd,Td] = schurAd(Ad,Q,T);
  Qtest = @(x) Q + x*Qd;
  Ttest = @(x) T + x*Td;
  Atest = @(x) Qtest(x)*Ttest(x)*Qtest(x)';
  normA = norm(A,'fro');
  ErrQ = @(x) Qtest(x)'*Qtest(x) - eye(20);
  ErrA = @(x) Atest(x) - (A + x*Ad);
  x = 10.^(-10:.1:0);
  ErrQRes = [];
  ErrARes = [];
  for i=1:length(x)
    ErrQRes = [ErrQRes, norm(ErrQ(x(i)),'fro')];
    ErrARes = [ErrARes, norm(ErrA(x(i)),'fro')./normA];
  end
  loglog(x,ErrQRes,'b;Q;',x,ErrARes,'r;A;')
end

#test the trainBust function
if(exist('tTrain','var') && tTrain == true)

  if(~exist('middle','var'))
    middle = 'b';
  end
  s = 8;
  A = hessen(randn(100));
  x = 100*randn(s,1);
  d = rand(s,1);
  ep = 1e-2;

  [~,~,sp1,~] = trainBust(A,x+ep*d,[],false,false,middle);
  [~,~,sp2,~] = trainBust(A,x-ep*d,[],false,false,middle);
  [~,deriv,sp0,~] = trainBust(A,x,d,false,false,middle);

  errR = (sp1 - sp0)/ep - deriv;
  errL = (sp0 - sp2)/ep - deriv;
  errC = .5*(sp1 - sp2)/ep - deriv;
  #The following vas inspired by noticing that with small ep (~1e-8) the
  #values in the spikes would swap.
  #Since the ordering of the spikes is arbitrary/can be sorted,
  #swapping 2 values should not really be represented by a discontinuity
  #Furthermore, in the formation of householder vectors, there are arbitrary (...you get my point...)
  #signs introduced, so a switch from plus to minus should not represent a discontinuity either - 
  #especially since the differentiation code assumes that the sign is a constant, not a function
  oddErrApp = .5*(sort(abs(sp1)) - sort(abs(sp2)))/ep;
  errOdd = sort(oddErrApp) - sort(deriv);

  printf('Ther errors were as follows:\nRHS (mean) = %f\nerrC (mean) = %f\nerrC (mean of log) = %f\nMax |errC| = %f\nsorted err max < %f\n',...
  norm(errR,1)/length(deriv), norm(errC,1)/length(deriv),...
  sum(log10(abs(errC) + eps/2))/length(errC), max(abs(errC)), max(abs(errOdd))...
  );

  #Justification, sorting locations is fine, all matrices are constant permutation matrices so doesn't effect
  #derivative except to permute it anyways.
  #The absolute value is justified because a switch in signs is not accounted for in the differentiation
  #(since arbitrary sign switches in householder formations are assumed to be constant)
 
end

if(exist('ltTrain','var') && ltTrain == true)
  s = 8;
  middle = 'm';
  
  Ex = 0;
  Ex2 = 0;
  V = [];
  D = [];
  O = [];
  for i=1:100
    A = hessen(randn(100));
    x = 100*randn(s,1);
    d = rand(s,1);
    d /= norm(d);
    ep = 1e-2;

    [~,~,sp1,~] = trainBust(A,x+ep*d,[],false,false,middle);
    [~,~,sp2,~] = trainBust(A,x-ep*d,[],false,false,middle);
    [~,deriv,sp0,~] = trainBust(A,x,d,false,false,middle);
    
    oddErrApp = .5*(sort(abs(sp1)) - sort(abs(sp2)))/ep;
    errOdd = log10(abs(sort(oddErrApp) - sort(deriv)));
    V = [V, errOdd];
    D = [D, sort(deriv)];
    O = [O, sort(oddErrApp)];
    Ex += sum(errOdd)/length(errOdd);
    Ex2 += norm(errOdd)^2/length(errOdd);
  end

  Ex /= 100;
  Ex2 /= 100;

  fprintf('Mean: %f\nStd: %f\n',Ex,sqrt(Ex2 - Ex^2));
  %hist(log10(abs((O - D)./(D+sqrt(eps)))))
  figure()
  hist(V)

end

#Tests analogous to tSchur2 which are included in adEffect.tex
if(exist('tradTrain','var') && tradTrain == true)
  #to generate the following without actually performing the
  #Schur decomposition step, go into trainBust and comment out the
  #line calling schurComp, then run this (figure 1a in adEffect)

  if(~exist('newA','var') || newA == true)
    A = hessen(randn(100)); 
    s = randn(10,1); 
    d = rand(10,1); 
    d /= norm(d);
  end
  [B,~,~,INFO] = trainBust(A,s,d);
  x = 10.^(-10:1/3:0);
  ers = zeros(length(x),1);
  for i=1:length(x); 
    [OUT{i},~,~,~] = trainBust(A,s+x(i)*d,[]);
    ers(i) = [norm(B + x(i)*INFO.dots{1}.H - OUT{i},'fro')]; 
  end
  figure();
  loglog(x,ers,'-;;','linewidth',4)
  set(gca(),'fontsize',16);
  xlabel('\epsilon','fontsize',20);
  #note: epsilon appears greek in windows, but prints to png as english word
  #print('notes/name.png','-dpng')

  #for figure 2, run this a couple of times to find satisfactory
  #results, turn on hold to make both appear on same plot
  #It took me 3 tries to get those two matrices. Other attemps yielded
  #results similar to one or the other of them.

end
