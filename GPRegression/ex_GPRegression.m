% (script)
% Example: GP regression 1D
%
% Use for lecture:
% 1) Show prior distribution only
% [comment 'return']
% 2) Show posterior (click through observations, iteratively add)
% [set flag_iterData = false]
% 3) vary hyperparameters
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
% [11.11.18, ST]    minor updates
% [27.01.17, ST]    minor update
% [24.01.17, ST]    plot prior/posterior next to each other, some added
%                   controls for plotting
% [22.01.17, ST]    minor updates
% [21.01.17, ST]    first version



%% Configuration
close all;
clear;

addpath([pwd,'/utils/']);

%%%% LUCAS RATH - BEGIN CHANGE
n_test = 100;
%%%% LUCAS RATH - END CHANGE

% Testpoints: these points will be used for plotting all functions

%%%% LUCAS RATH - BEGIN CHANGE
% x_plot = linspace(-4,4)';  % test points (plot functions at these points)
x_plot = linspace(-4,4,n_test);  % <D,n> test points (plot functions at these points)
% n_test = length(x_plot);
% n_test = size(x_plot,2);
%%%% LUCAS RATH - END CHANGE

% plotting
cmp = get(groot,'DefaultAxesColorOrder'); % new colormap
c_gray = 0.4*ones(1,3);

lineWidth1 = 1.;
lineWidth2 = 2;


pos_fig = [50 50 1024 420];
%pos_fig = [50 50 605 540];  % one graph only

% y-lim function
ylim_f = [-3,3];


%% User control
% What to show:
% 1) ell: 0.2, 2 (best: 2)
% 2) sigma_f: 0.05 (best: 1)
% 3) sigma_n: 0.2 (best: 0.05)

% Training points (order matters for plotting)
%X = -4:1:4;
X = [1 2 3 -1 -3 -4 4 0];

% Flag: iteratively add data?
flag_iterData = false;

% Show samples or not?
flag_noSamples = false;

% Hyperparameters
sigma_f = 1;      % output variance (std); default/best: 1
ell = 2;            % length scale; default: 1, best: 2
sigma_n = 0.05;     % STD of measurement noise; default/best: 0.05



%% Get data
% Training data
n_train = length(X);
Y = truthFcn(X)' + sigma_n*randn(n_train,1);

% plot truth function and data
figure;
set(gcf,'Position',pos_fig);
plot(x_plot,truthFcn(x_plot),'color',cmp(1,:));
set(gca,'ylim',ylim_f);
grid on;
hold on;
plot(X,Y,'.r','markersize',15);
drawnow;


%% Prior
% Mean function
%gp.mu = @(a) 0;    % one data point

%%%% LUCAS RATH - BEGIN CHANGE
% gp.mu = @(a) (zeros(size(a,1),1));  % multiple data points
gp.mu = @(a) (zeros(size(a,2),1));  % <n,1>  multiple data points
% gp.mu = @(a) (linspace(1,5,size(a,2)))';  % <n,1> multiple data points

% Y = [X;ones(1,n_train)]' * W
% least-square-fit =>  W = [X;ones(1,n_train)]'\Y
% wprior = [X;ones(1,size(X,2))]'\Y;
% gp.mu = @(a) [a;ones(1,size(a,2))]' * wprior;
%%%% LUCAS RATH - END CHANGE

% Kernel: squared exponential (for scalar input)
%gp.k = @(a,b) sigma_f^2 * exp(-(a-b)^2/(2*ell^2));    % one data point
gp.K = @(a,b) sigma_f^2 * exp(-bsxfun(@minus,a,b').^2/(2*ell^2));    % <na,nb> multiple data points


% GRV prior at test points
x = x_plot;
m_x = gp.mu(x);
K_xx = gp.K(x,x);

% Plot
figure;
subplot(1,2,1);
set(gcf,'Position',pos_fig);
[h_fig] = MyGPplot(gcf,m_x,K_xx,x,[],[],[],flag_noSamples);
plot(x_plot,truthFcn(x_plot),'color',cmp(1,:),'linewidth',lineWidth1); % truth function

grid on;
xlabel('x');
ylabel('f(x)');
set(gca,'ylim',ylim_f);
title('Prior GP');

%return;


%% Posterior

%figure;
subplot(1,2,2);
set(gcf,'Position',pos_fig);
title('Posterior GP');

% Add one training point after the other?
if flag_iterData
    % iteratively add data
    iters = 1:length(X);
else
    % don't
    iters = length(X);
end
for n_sel=iters
    % data points considered
    X_sel = X(1:n_sel);
    Y_sel = Y(1:n_sel);
    
    % NOTE: this is not the most numerically effective implementation (esp.
    % matrix inversion should not be done this way), but meant to closely 
    % resemble the lecture.  See, e.g. book by Rasmussen/Williams, p. 19
    % for a numerically more effective implementation.
    m_X = gp.mu(X_sel);
    k_XX = gp.K(X_sel',X_sel');
    aux1 = inv(k_XX + sigma_n^2*eye(n_sel));

    %%%% LUCAS RATH - BEGIN CHANGE
    % k_xX = gp.K(x,X_sel);
    k_xX = gp.K(x',X_sel');
    %%%% LUCAS RATH - END CHANGE
    
    mean_post = m_x + k_xX*aux1*(Y_sel-m_X);
    K_post = K_xx - k_xX*aux1*k_xX';
    
    % plot
    cla;
    [h_fig] = MyGPplot(gcf,mean_post,K_post,x,[],[],[],flag_noSamples);
    %
    % Truth function
    plot(x_plot,truthFcn(x_plot),'Color',cmp(1,:),'linewidth',lineWidth1); % truth function
    %
    % Add measurements with error bars
    errorbar(X_sel,Y_sel,2*sigma_n*ones(size(Y_sel)),'o','Color',cmp(1,:),'Linewidth',lineWidth1);
    set(gca,'ylim',ylim_f);
    xlabel('x');
    ylabel('f(x)');

    if flag_iterData
        pause;
    end
end

return;

