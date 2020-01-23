clear all; close all; clc;

var_w   = 1e-8;
var_f   = 1e4;              % output variance
M       = diag(1e0);     % length scale
var_n   = var_w;           % measurement noise variance
maxsize = 100;             % maximum number of points in the dictionary

% create GP object
gp = GP(1, 1, var_f, var_n, M, maxsize);

%% true function
Xt = [-4 -3 -2 -1 0 1 2 3 4];
Yt = [-1 0.5 1.5 1.5 1 0.8 1 1.5 2]';
pp = spline(Xt,Yt');
truefun = @(x) ppval(pp, x)';

%% measure data
X = Xt([1 2 3 5 7 8]);
Y = truefun(X) + sqrt(var_n)*randn(length(X),1);

gp.add(X,Y);


%% plot
close all;

% test 1
gp.M = 1e-2;
gp.var_f = 2;
gp.var_n = var_w;
gp.updateModel();
gp.plot1d(truefun)
xlim([-6 5]); ylim([-4 4]);

% test 2
gp.M = 50;
gp.var_f = 100;
gp.var_n = var_w;
gp.updateModel();
gp.plot1d(truefun)
xlim([-6 5]); ylim([-4 4]);

% test 3
gp.M = 1e-3;
gp.var_f = 1e-2;
gp.var_n = var_w;
gp.updateModel();
gp.plot1d(truefun)
xlim([-6 5]); ylim([-4 4]);

% test 4
gp.M = 50;
gp.var_f = 1e7;
gp.var_n = var_w;
gp.updateModel();
gp.plot1d(truefun)
xlim([-6 5]); ylim([-4 4]);

% test 4
gp.M = 100^2;
gp.var_f = 100^2;
gp.var_n = 1^2;
gp.updateModel();
gp.plot1d(truefun)
xlim([-6 5]); ylim([-4 4]);

%% optimize and plot

gp.M = 0.1^2;
gp.var_f = 0.1^2;
gp.var_n = 1e-8;

% gp.optimizeHyperParams('ga');
gp.optimizeHyperParams('fmincon');
gp.plot1d(truefun)
xlim([-6 5]); ylim([-4 4]);
