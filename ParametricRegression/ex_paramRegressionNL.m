% (script)
% Example: parametric regression (NONLINEAR model) 1D
%
% Use for lecture:
% 1) Show posterior for linear feature (first with data 2:4, then -4:4)
% [change features]
% 2) Show posterior for affine
% etc.
%
% --
% Matlab tutorial on Gaussian process regression.
%
% Copyright 2017-19 Max Planck Society. All rights reserved.
% 
% Sebastian Trimpe
% Max Planck Institute for Intelligent Systems
% trimpe(at)is.mpg.de
%
% Revision history
% [22.01.17, ST]    minor updates
% [21.01.17, ST]    first version
% [11.11.18, ST]    minor updates


%% Configuration
close all;
clear;

% plotting
cmp = get(groot,'DefaultAxesColorOrder'); % new colormap
c_gray = 0.4*ones(1,3);

lineWidth1 = 1.5;
lineWidth2 = 2;

pos_fig = [50 50 605 540];
%pos_fig = [50 556 605 540];

% y-lim function
ylim_f = [-1.5, 3];


%% User control
% What to show:
% 1) linear + -4:1:4
% 2) affine
% 3) quadratic
% 4) x.^3

% Features
% phi = @(x) x;   % linear
%phi = @(x) [ones(1,length(x)); x];   % affine
%phi = @(x) [ones(1,length(x)); x; x.^2];   % quadradic
%phi = @(x) [ones(1,length(x)); x; x.^2; x.^3];
phi = @(x) [ones(1,length(x)); x; x.^2; x.^3; x.^4; x.^5];

% Training points
X = -4:1:4;
% X = [2:4];


%% Get data
% Training data
sigma_n = 0.1;      % STD of measurement noise

n_train = length(X);
Y = truthFcn(X)' + sigma_n*randn(n_train,1);

% plot truth function and data
figure;
set(gcf,'Position',pos_fig);
x_plot = linspace(-4,4)';
plot(x_plot,truthFcn(x_plot));
set(gca,'ylim',ylim_f);
hold on;
grid on;
plot(X,Y,'.r','markersize',15);
%drawnow;


%% Prior
N = length(phi(0));     % number of features

w_mean = zeros(N,1);
w_var = eye(N);
w_std = chol(w_var)';


%% Posterior

% posterior 
phi_X = phi(X);
w_mean_post = w_mean + ...
    w_var*phi_X*inv(sigma_n^2*eye(n_train)+phi_X'*w_var*phi_X)*(Y-phi_X'*w_mean);
w_var_post = w_var - w_var*phi_X*inv(sigma_n^2*eye(n_train)+phi_X'*w_var*phi_X)*phi_X'*w_var;


% Plot 
figure;
set(gcf,'Position',pos_fig);

% Use GP-plotting function
xx = linspace(-4,4);
%xx = linspace(-6,6);
phi_xx = phi(xx);
m = phi_xx'*w_mean_post;
V = phi_xx'*w_var_post*phi_xx;
[h_fig] = MyGPplot(gcf,m,V,xx,[],[],[],true);
grid on;
xlabel('x');
ylabel('f(x)');

% Add data
hold on;
plot(X,Y,'.r','markersize',12');

set(gca,'ylim',ylim_f);
