function [y]=truthFcn(x)
% True function to be approximated.
%
% --
% Matlab tutorial on Gaussian process regression.
%
% Copyright 2017 Max Planck Society. All rights reserved.
% 
% Sebastian Trimpe
% Max Planck Institute for Intelligent Systems
% Autonomous Motion Department
% strimpe(at)tuebingen.mpg.de
%
% Revision history
% [24.06.16, ST]    first version
% [23.01.17, ST]    minor update


%% Determine coefficients of polynomial.
% Truth function is a spline determined by the following points
X = [-4 -3  -2  -1  0 1   2 3   4];
Y = [-1 0.5 1.5 1.5 1 0.8 1 1.5 2];

pp = spline(X,Y);


%% Evaluate function
y = ppval(pp, x);