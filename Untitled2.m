close all;
clc;
clear;

% register the generator
setpaths;

m = gbs_Matrix('m', 2, 4); % three 2D measurements
n = gbs_Vector('n', 7);
d = gbs_Vector('d', 3); % three known 3D distances
syms a1 a2 a3; % three unknown lambda factors
eq(2) = a1^2+a2^2-m(1,1)*a1*a2-m(1,2)*a1+m(1,3)*a2+m(1,4)-d(1); % three polynomial equations
eq(3) = a1^2+a3^2-m(2,1)*a1*a3-m(2,2)*a1+m(2,3)*a3+m(2,4)-d(2);
eq(1) = a1^2-n(1)*a1*a2-n(2)*a1*a3+n(3)*a2*a3+n(4)*a1+n(5)*a2+n(6)*a3+n(7)-d(3); 
unknown={'a1','a2','a3'};
known={'m11','m12','m13','m14','m21','m22','m23','m24','n1','n2','n3','n4','n5','n6','n7','d1','d2','d3'};
[res,export]=gbs_CreateCode('NSP3P_solver',eq,known,unknown);