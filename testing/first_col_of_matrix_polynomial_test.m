
% Show off an efficient way to generate the first column of the matrix
% polynomial

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


shifts = [2, 3, 5, 7.8, 9, 17, 25];



% Generate first column of (A - s1) (A - s2) (A - s3)...
%
% Note for C impl: use transposed M so that v can be directly strided into M
% (saves a copy, and the end+1'th element is allocated for free every time).

ns = length(shifts);
Asmall = A(1:ns, 1:ns);
straggler = A(ns+1, ns);

% compute all but last element using matrix-matrix multiplies
M = eye(ns);
for s = shifts
	Asmall = Asmall - s*eye(ns);
	M = M * Asmall;
	Asmall = Asmall + s*eye(ns);
end
v = M(:, 1);

% compute last element (product of subdiagonal entries)
v(end+1) = straggler * prod(diag(Asmall, -1));



% The old-fashioned way

H = eye(size(A));

for s = shifts
	H = H * (A - eye(size(H)) * s);
end

err = v - H(1:length(v), 1);


v
sum_square_error = err' * err

