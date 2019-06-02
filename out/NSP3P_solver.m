% Generated using GBSolver generator Copyright Martin Bujnak,
% Zuzana Kukelova, Tomas Pajdla CTU Prague 2008.
% 
% Please refer to the following paper, when using this code :
%     Kukelova Z., Bujnak M., Pajdla T., Automatic Generator of Minimal Problem Solvers,
%     ECCV 2008, Marseille, France, October 12-18, 2008

function [a1 a2 a3] = NSP3P_solver(m11, m12, m13, m14, m21, m22, m23, m24, n1, n2, n3, n4, n5, n6, n7, d1, d2, d3)

	% precalculate polynomial equations coefficients
	c(1) = 1;
	c(2) = -n1;
	c(3) = -n2;
	c(4) = n3;
	c(5) = n4;
	c(6) = n5;
	c(7) = n6;
	c(8) = n7-d3;
	c(9) = 1;
	c(10) = -m11;
	c(11) = 1;
	c(12) = -m12;
	c(13) = m13;
	c(14) = m14-d1;
	c(15) = 1;
	c(16) = -m21;
	c(17) = 1;
	c(18) = -m22;
	c(19) = m23;
	c(20) = m24-d2;

	M = zeros(24, 32);
	ci = [17, 40, 87, 110, 181, 294, 317, 388, 529];
	M(ci) = c(1);

	ci = [41, 64, 111, 134, 205, 318, 341, 412, 553];
	M(ci) = c(2);

	ci = [113, 136, 183, 206, 253, 390, 413, 460, 601];
	M(ci) = c(3);

	ci = [137, 160, 207, 230, 277, 414, 437, 484, 625];
	M(ci) = c(4);

	ci = [329, 352, 399, 422, 469, 534, 557, 604, 673];
	M(ci) = c(5);

	ci = [353, 376, 423, 446, 493, 558, 581, 628, 697];
	M(ci) = c(6);

	ci = [425, 448, 471, 494, 517, 606, 629, 652, 721];
	M(ci) = c(7);

	ci = [569, 592, 615, 638, 661, 678, 701, 724, 745];
	M(ci) = c(8);

	ci = [21, 92, 115, 186, 297, 320, 391, 530];
	M(ci) = c(9);

	ci = [45, 116, 139, 210, 321, 344, 415, 554];
	M(ci) = c(10);

	ci = [69, 140, 163, 234, 345, 368, 439, 578];
	M(ci) = c(11);

	ci = [333, 404, 427, 474, 537, 560, 607, 674];
	M(ci) = c(12);

	ci = [357, 428, 451, 498, 561, 584, 631, 698];
	M(ci) = c(13);

	ci = [573, 620, 643, 666, 681, 704, 727, 746];
	M(ci) = c(14);

	ci = [48, 95, 118, 300, 323, 394, 531];
	M(ci) = c(15);

	ci = [144, 191, 214, 396, 419, 466, 603];
	M(ci) = c(16);

	ci = [240, 263, 286, 468, 491, 514, 651];
	M(ci) = c(17);

	ci = [360, 407, 430, 540, 563, 610, 675];
	M(ci) = c(18);

	ci = [456, 479, 502, 612, 635, 658, 723];
	M(ci) = c(19);

	ci = [600, 623, 646, 684, 707, 730, 747];
	M(ci) = c(20);


	Mr = rref(M);  % replace me with a MEX

	A = zeros(8);
	amcols = [32 31 30 29 28 27 26 22];
	A(1, 4) = 1;
	A(2, 7) = 1;
	A(3, :) = -Mr(23, amcols);
	A(4, :) = -Mr(22, amcols);
	A(5, :) = -Mr(20, amcols);
	A(6, :) = -Mr(18, amcols);
	A(7, :) = -Mr(17, amcols);
	A(8, :) = -Mr(11, amcols);

	[V D] = eig(A);
	sol =  V([4, 3, 2],:)./(ones(3, 1)*V(1,:));

	if (find(isnan(sol(:))) > 0)
		
		a1 = [];
		a2 = [];
		a3 = [];
	else
		
		I = find(not(imag( sol(1,:) )));
		a1 = sol(1,I);
		a2 = sol(2,I);
		a3 = sol(3,I);
	end
end
