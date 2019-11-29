%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -
%------------------------------------------------------------------

clear vars; close all; clc;

% GP hyperparameters
sigmaf  = 10;          % output variance (std)
lambda  = eye(2)*10^2; % length scale
sigman  = 1;           % STD of measurement noise
maxsize = 100;         % maximum number of points in the dictionary

% create GP object
gp = GP(sigmaf, sigman, lambda, maxsize);

% true function
truthfun = @(x) x'*[1;2] - 20;

% generate sampled data
nsamples = 20;
xdata = rand(2,nsamples)*10 + 10;
ydata = truthfun(xdata) + sigman*randn(size(xdata,2),1);

% add sampled data to gp dictionary
gp.add(xdata,ydata)
% gp.add(xdata,[ydata,ydata])   % test 2D output case

% evaluate at some arbitrary points to check if GP works
x = [25:35; 25:35];
[muy, covary] = gp.eval(x);

% plot prediction bias and variance
gp.plot2d(truthfun, [-5,40],[-5,40])
