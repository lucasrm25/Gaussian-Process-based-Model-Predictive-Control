% (script)
% Example: parametric regression (linear model) 1D
%
% Use for lecture:
% 1) Show prior distribution only (click through samples)
% [comment 'return']
% 2) Show posterior (click through samples)
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
% [20.01.17, ST]    first version
% [22.01.17, ST]    minor updates
% [11.11.18, ST]    minor updates


%% Configuration
close all;
clear;

addpath([pwd,'/utils/']);


% plotting
cmp = get(groot,'DefaultAxesColorOrder'); % new colormap
c_gray = 0.4*ones(1,3);

lineWidth1 = 1.5;
lineWidth2 = 2;

pos_fig = [50 50 1024 420];
%pos_fig = [50 676 1024 420];

% y-lim function
ylim_f = [-8, 8];



%% Get data
% Training data
sigma_n = 1;      % STD of measurement noise
data_X = [2 3]';
n_train = length(data_X);
%Y = truthFcn(X) + sigma_n*randn(n_train,1);
data_Y = data_X*0.5;  % for simplicity: hard-code for function f = 0.5*x

% % plot truth function and data
% x_plot = linspace(-4,4)';
% plot(x_plot,truthFcn(x_plot));
% hold on;
% grid on;
% plot(data_X,data_Y,'.r','markersize',12);


%% Prior
w_mean = 0;
w_var = 1;
w_std = sqrt(w_var);

% Plot PDF of w
figure;
set(gcf,'Position',pos_fig);
[ww,pdf_w] = plot_paramEx(w_mean,w_std,[],5,ylim_f);
legend('prior');

%return;


%% Posterior
% measurement
Y = data_Y(1);
X = data_X(1);

% likelihood
w_lik = normpdf(Y, X*ww, sigma_n);

% posterior 
w_mean_post = w_mean + w_var*X*inv(sigma_n^2+X*w_var*X)*(Y-X*w_mean);
w_var_post = w_var - w_var*X*inv(sigma_n^2+X*w_var*X)*X*w_var;
% Checked: gives same result:
% pdf_w_post = pdf_w.*w_lik;
% pdf_w_post_normal = trapz(ww,pdf_w_post);
% pdf_w_post = pdf_w_post/pdf_w_post_normal;


% Plot 
figure;
set(gcf,'Position',pos_fig);

subplot(1,2,1);
plot(ww,pdf_w,'color',cmp(2,:),'linewidth',lineWidth1); hold on;
plot(ww,w_lik,'color',cmp(3,:),'linewidth',lineWidth1);
%plot(ww,normpdf(ww,w_mean_post,sqrt(w_var_post)),'color',cmp(2,:));

data.x = X;
data.y = Y;
[~,pdf_w_post] = plot_paramEx(w_mean_post,sqrt(w_var_post),data,5,ylim_f);

subplot(1,2,1);
legend('prior','likelihood','posterior');


%% Posterior (2nd measurement)
% measurement
Y = data_Y(2);
X = data_X(2);

% likelihood
w_lik = normpdf(Y, X*ww, sigma_n);

% posterior 
w_mean_post2 = w_mean_post + w_var*X*inv(sigma_n^2+X*w_var*X)*(Y-X*w_mean_post);
w_var_post2 = w_var_post - w_var_post*X*inv(sigma_n^2+X*w_var_post*X)*X*w_var_post;


% Plot 
figure;
set(gcf,'Position',pos_fig);

subplot(1,2,1);
plot(ww,pdf_w_post,'color',cmp(2,:),'linewidth',lineWidth1); hold on;
plot(ww,w_lik,'color',cmp(3,:),'linewidth',lineWidth1);
%plot(ww,normpdf(ww,w_mean_post,sqrt(w_var_post)),'color',cmp(2,:));

data.x = X;
data.y = Y;
plot_paramEx(w_mean_post2,sqrt(w_var_post2),data,5,ylim_f);

subplot(1,2,1);
legend('prior','likelihood','posterior');
