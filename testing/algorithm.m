
addpath ../src

N=13;
A=[2.81014,-6.48911,-15.18678,21.77452,5.07491,-33.61697,7.69500,-3.21216,41.49093,16.75451,37.52472,-4.91861,30.00529;
-79.66993,32.07671,-37.41973,27.72298,3.11131,-59.35638,-17.86218,18.53775,72.09095,3.57814,158.46592,284.03035,333.00127;
0,17.26496,33.19817,-38.30290,-25.19686,31.03556,-58.71160,35.69149,-70.85779,-69.20413,73.02552,-81.81205,-153.14354;
0,0,-4.70294,19.82071,17.59862,-9.54619,22.05326,-23.91293,23.02401,50.87808,-87.43072,41.81270,130.30457;
0,0,0,6.46713,17.41596,-7.30266,-1.34361,0.41457,3.50943,-5.91071,34.33494,15.76940,8.73907;
0,0,0,0,15.48683,20.92098,26.03202,-17.96809,7.73529,46.79353,-112.27704,-9.34287,71.89896;
0,0,0,0,0,11.00435,26.79193,4.71123,-4.04375,-14.29687,16.39901,18.17265,-39.43072;
0,0,0,0,0,0,5.19289,19.11160,3.20148,16.21861,-28.96126,-0.84516,29.94569;
0,0,0,0,0,0,0,7.72795,10.10828,9.25127,-2.09658,-7.79144,3.98539;
0,0,0,0,0,0,0,0,3.19594,12.23062,-9.39841,4.20665,6.76808;
0,0,0,0,0,0,0,0,0,-6.72820,22.71813,-19.60648,-28.18783;
0,0,0,0,0,0,0,0,0,0,-1.11828,25.13190,9.08088;
0,0,0,0,0,0,0,0,0,0,0,-20.40095,-4.33512];



%N=19;
%r = rand([N, N]);
%A = hess(r*diag(list_primes(N))*inv(r));



% Generate first column of (A - s1) (A - s2) (A - s3)...
%
% Note for C impl: use transposed M so that v can be directly strided into M
% (saves a copy, and the end+1'th element is allocated for free every time!).
function [v] = firstColMatrixPolynomial(Asmall, last_subdiagonal, shifts)
	ns = length(shifts);

	% compute all but last element using matrix-matrix multiplies
	M = eye(ns);
	for s = shifts
		Asmall = Asmall - s*eye(ns);
		M = M * Asmall;
		Asmall = Asmall + s*eye(ns);
	end
	v = M(:, 1);

	% compute last element (product of first ns subdiagonal entries)
	v(end+1) = last_subdiagonal * prod(diag(Asmall, -1));
end



% Generate householder reflector which annihilates
%
%     [. . . . . . . .]^T
%
% into
%
%     [. 0 0 0 0 0 0 0]^T
%
function [Q] = householderAnnihilator(v)
	v(1) += sgn(v(1)) * sqrt(v'*v);
	v = v / norm(v);
	Q = eye(length(v)) - 2 * v * v';
end
function y = sgn(x)
	y = sign(x);

	if (x == 0)
		y = 1;
	end
end



% submatrix product and "trailing update", i.e. multiply
%
% [ I  0  0 ]       [ I  0   0 ]   [    A11    A12 Q'    A13 ]
% [ 0  Q  0 ] * A * [ 0  Q'  0 ] = [  Q A21  Q A22 Q'  Q A23 ]
% [ 0  0  I ]       [ 0   0  I ]   [    A31    A32 Q'    A33 ]
%
% where the number of rows (and cols) in A11 is offset. This (potentially0
% allows for parallelization later on
%
% TODO - specialize this for Hessenberg matrices? can shave some multiplies
function [M] = submatrixProduct(M, offset, Q)
	idx = (offset + 1):(offset + length(Q));
	M(idx, :) = Q * M(idx, :);
	M(:, idx) = M(:, idx) * Q';
end



% backwards matrix reference / reference from Bottom Right
function [M] = brref(A,r,c)
	l = length(A);
	M = rot90(rot90(A(l+1-r,l+1-c)));
end



% transpose about the antidiagonal ("wrong" diagonal)
function [M] = adtranspose(M)
	M = rot90(rot90(transpose(M)));
end











shifts = [41 37 31 29]
ns = length(shifts);

%% create bulge on top
%A = submatrixProduct(A, 0, ...
%      householderAnnihilator( ...
%        firstColMatrixPolynomial(A(1:ns, 1:ns), A(ns+1, ns), shifts)   ));

% create bulge on bottom
A = submatrixProduct(A, N-ns-1, ...
      householderAnnihilator( ...
        firstColMatrixPolynomial(A(1:ns, 1:ns), A(ns+1, ns), fliplr(shifts))   ));

% chase bottom bulge
for i = 1:(N-ns-1)
	j=i+1;
	A = submatrixProduct(A, N-ns-j, ...
		  adtranspose(householderAnnihilator(adtranspose(brref(A, i, j:(j+ns)))))   );
end

disp "after chase"
eig(A)'

% spikes
r = 1:ns;
[Q, ~] = schur(A(r, r), "real");
A = submatrixProduct(A, 0, Q');

% YIPPIE
disp("eig of top left 4x4, bottom right 9x9");
r = 1:ns;     eig(A(r,r))'
r = (ns+1):N; eig(A(r,r))'
A = hess(A)

disp("eig of top left 4x4, bottom right 9x9");
r = 1:ns;     eig(A(r,r))'
r = (ns+1):N; eig(A(r,r))'

