close all;
clc;
clear;

% register the generator
setpaths;

m = gbs_Matrix('m', 4, 6); % three 2D measurements

syms a b c; % three unknown lambda factors
s=[a^2,a,b,b*c^2,b*c,1].';
eq=m*s;
unknown={'a','b','c'};
known={'m11','m12','m13','m14','m15','m16',...
       'm21','m22','m23','m24','m25','m26',...
       'm31','m32','m33','m34','m35','m36',...
       'm41','m42','m43','m44','m45','m46'};
[res,export]=gbs_CreateCode('solver',eq,known,unknown);