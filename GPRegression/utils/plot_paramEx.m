function [ww,pdf_w] = plot_paramEx(w_mean,w_std,data,n_sample,my_ylim)
% Plotting function for parametric example.
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
% [20.01.17, ST]    first version

if nargin<3
    data = [];
end;
if nargin<4
    n_sample = 0;
end;
if nargin<5
    my_ylim = [];
end;

%% Config
cmp = get(groot,'DefaultAxesColorOrder'); % new colormap
c_gray = 0.4*ones(1,3);

lineWidth1 = 1.5;
lineWidth2 = 2;

% Pause when plotting samples
flag_pause = true;


%% Plot PDF of w
subplot(1,2,1);
ww = linspace(-4,4,1000);
pdf_w = normpdf(ww,w_mean,w_std);
plot(ww,pdf_w,'color',cmp(1,:),'linewidth',lineWidth1); hold on;
grid on;
xlabel('w');
ylabel('p(w)');
title('Distribution of weights');


%% Plot function f(x)
subplot(1,2,2);

% Use GP-plotting function
xx = linspace(-4,4);
m = (w_mean * xx)';
V = xx'*w_std^2*xx;
[h_fig] = MyGPplot(gcf,m,V,xx,[],[],[],true);
grid on;
xlabel('x');
ylabel('f(x)');
title('Distribution of functions');

if ~isempty(my_ylim)
    set(gca,'ylim',my_ylim);
end;
if ~isempty(data)
    plot(data.x,data.y,'.r','markersize',12');
end;


%% Add samples
% samples
%n_sample = 10;   % number of samples
w_s = w_mean + w_std*randn(n_sample,1);
for i=1:n_sample
    if flag_pause
        pause;
    end;
    subplot(1,2,1);
    plot([w_s(i), w_s(i)], [0 normpdf(w_s(i),w_mean,w_std)],'--','color',c_gray);
    subplot(1,2,2);
    plot(xx,w_s(i)*xx,'--','color',c_gray);
end;

