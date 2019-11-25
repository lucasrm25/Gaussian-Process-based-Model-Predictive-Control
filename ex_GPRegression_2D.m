clear all; close all; clc;

% GP hyperparameters
sigmaf  = 10;        % output variance (std)
lambda  = 10;        % length scale
sigman  = 0.05;     % STD of measurement noise
maxsize = 100;      % maximum number of points in the dictionary

% create GP object
gp = GP(sigmaf, sigman, lambda, maxsize);

% generate sampled data
truthfun = @(x) x'*[1;2] - 20;
xdata = rand(2,10)*10 + 10;
ydata = truthfun(xdata) + sigman*randn(size(xdata,2),1);

% add sampled data to gp dictionary
gp.add(xdata,ydata)

% evaluate at some arbitrary points to check if algo. works
x = [25:35; 25:35];
[muy, covary] = gp.eval(x);

% plot prediction bias and variance
gp.plot2d(-5:40,-5:40, truthfun)