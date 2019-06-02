close all;
clc;
clear;
% register the generator
setpaths;

g = gbs_Vector('g', 4);
syms x y; % three unknown lambda factors
eq(1) = x^2*g(1)-y^2*g(1)+x*y*g(2)+x*g(3)-y*g(4);
eq(2) = x^2+y^2-1;
unknown={'x','y'};
known={'g1','g2','g3','g4'};
[res,export]=gbs_CreateCode('solver',eq,known,unknown);