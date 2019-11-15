% (script)
% Plot data
%
% --
% Matlab tutorial on Gaussian process optimization.
%
% Copyright 2017 Max Planck Society. All rights reserved.
% 
% Sebastian Trimpe
% Max Planck Institute for Intelligent Systems
% Autonomous Motion Department
% strimpe(at)tuebingen.mpg.de
%
% Revision history
% [20.01.17, ST]    first version


%% Configuration
close all;
clear;



%% Get data
% Training data
sigma_n = 0.05;      % STD of measurement noise
%X = [-4 -3 -2 -1 1 2 3]';
X = (-4:1:4)';
n_train = length(X);
Y = truthFcn(X) + sigma_n*randn(n_train,1);

% plot truth function and data
x_plot = linspace(-4,4)';
plot(x_plot,truthFcn(x_plot));
hold on;
grid on;
plot(X,Y,'.r','markersize',12);
