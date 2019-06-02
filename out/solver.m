% Generated using GBSolver generator Copyright Martin Bujnak,
% Zuzana Kukelova, Tomas Pajdla CTU Prague 2008.
% 
% Please refer to the following paper, when using this code :
%     Kukelova Z., Bujnak M., Pajdla T., Automatic Generator of Minimal Problem Solvers,
%     ECCV 2008, Marseille, France, October 12-18, 2008

function [x y] = solver(g1, g2, g3, g4)

	% precalculate polynomial equations coefficients
	c(1) = g1;
	c(2) = g2;
	c(3) = -g1;
	c(4) = g3;
	c(5) = -g4;
	c(6) = 0;
	c(7) = 1;
	c(8) = 1;
	c(9) = -1;

	M = zeros(6, 10);
	ci = [4, 9, 25];
	M(ci) = c(1);

	ci = [10, 15, 31];
	M(ci) = c(2);

	ci = [16, 21, 37];
	M(ci) = c(3);

	ci = [28, 33, 43];
	M(ci) = c(4);

	ci = [34, 39, 49];
	M(ci) = c(5);

	ci = [46, 51, 55];
	M(ci) = c(6);

	ci = [6, 11, 26];
	M(ci) = c(7);

	ci = [18, 23, 38];
	M(ci) = c(8);

	ci = [48, 53, 56];
	M(ci) = c(9);


	Mr = rref(M);  % replace me with a MEX

	A = zeros(4);
	amcols = [10 9 8 7];
	A(1, 3) = 1;
	A(2, :) = -Mr(6, amcols);
	A(3, :) = -Mr(5, amcols);
	A(4, :) = -Mr(3, amcols);

	[V D] = eig(A);
	sol =  V([3, 2],:)./(ones(2, 1)*V(1,:));

	if (find(isnan(sol(:))) > 0)
		
		x = [];
		y = [];
	else
		
		I = find(not(imag( sol(1,:) )));
		x = sol(1,I);
		y = sol(2,I);
	end
end
