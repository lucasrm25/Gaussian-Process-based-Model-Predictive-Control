function [h_fig] = MyGPplot(h_fig,m,V,z,colChoice,paramThrsh,flag_noShading,flag_noSamples)
% GP plot.  Largely based on Philipp Hennig's code, lecture etc., see [1].
%
%   h_fig   figure handle
%   m       GP mean
%   V       GP variance
%   z       grid points
%   colChoice   'r' for red, 'g' for green
% 
% [1] Philipp Hennig, "Animating Samples from Gaussian Distributions," 
% Technical Report No. 8 of the Max Planck Institute for Intelligent 
% Systems, September 2013.
%
% --
% Matlab tutorial on Gaussian process regression.
% Largely based on Philipp Hennig's code, lecture etc.
%
% Copyright 2016/17 Max Planck Society. All rights reserved.
% 
% Sebastian Trimpe
% Max Planck Institute for Intelligent Systems
% Autonomous Motion Department
% strimpe(at)tuebingen.mpg.de
%
% Revision history
% [??, ST]    first version
% [23.01.17, ST]    minor update

if nargin<5 || isempty(colChoice)
    colChoice = 'r';
end;
if nargin<6 || isempty(paramThrsh)
    paramThrsh = 0.4;   % "the smaller, the brighter", Philipp: 0.6
end;
if nargin<7 || isempty(flag_noShading)
    flag_noShading = false;
end;
if nargin<8 || isempty(flag_noSamples)
    flag_noSamples = false;
end;

%% Plot config
% Number of samples to be plot
N_sample = 3;

% Colors
dgr = [0,0.4717,0.4604]; % color [0,125,122]
dre = [0.4906,0,0]; % color [130,0,0]
lightdgr = [1,1,1] - 0.5 * ([1,1,1] - dgr);
lightdre = [1,1,1] - 0.5 * ([1,1,1] - dre);
% dgr2white =
% bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,0.6,2024)').^0.5,[1,1,1]-dgr));
% % 2024 doesn't work!!!
% dre2white = bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,0.6,2024)').^0.5,[1,1,1]-dre));
dgr2white = bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,paramThrsh,256)').^0.5,[1,1,1]-dgr));
dre2white = bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,paramThrsh,256)').^0.5,[1,1,1]-dre));

if colChoice=='g'
    col_d = dgr;
    col_l = lightdgr;
    col_2white = dgr2white;
elseif colChoice=='r'
    col_d = dre;
    col_l = lightdre;
    col_2white = dre2white;
else
    error('Color not specified.');
end;

% Lines
lineWidthMean = 2;
lineWidthStd = 0.8;
lineWidthSamples = 0.8;

% Number of std
n_std = 2;  % how many std to plot
multSTD = 1.5*n_std;    % 4,3;  2, 2.5; Multiply STD to get shading y-range



%% Calcs
N_test = length(z);
std_V = sqrt(diag(V));
% max_std_V_shading = max(std_V*multSTD+m);
% min_std_V_shading = min(-std_V*multSTD+m);
max_std_V_shading = max([std_V*multSTD+m; -std_V*multSTD+m]);
min_std_V_shading = min([std_V*multSTD+m; -std_V*multSTD+m]);
max_std_V = max([std_V*n_std+m; -std_V*n_std+m]);
min_std_V = min([std_V*n_std+m; -std_V*n_std+m]);

% Gaussian density
GaussDensity = @(y,m,v)(bsxfun(@rdivide,exp(-0.5*...
    bsxfun(@rdivide,bsxfun(@minus,y,m').^2,v'))./sqrt(2*pi),sqrt(v')));

% Samples from prior
s_prior = bsxfun(@plus, m, chol(V + 1.0e-8 * eye(N_test))' * randn(N_test,N_sample));

% Shading
%Y = linspace(-multSTD*max_std_V,multSTD*max_std_V,250)';
Y = linspace(min_std_V_shading,max_std_V_shading,250)';
zGD = linspace(z(1),z(end),300); % x-resolution of the density
P = GaussDensity(Y,interp1(z,m,zGD)',interp1(z,diag(V+eps),zGD)'); 


%% Plot
if ~isempty(h_fig)
    figure(h_fig);
end;

% Set colormap
colormap(col_2white);
%clf;
hold on; 
box on;
grid on;

% Shading
if ~flag_noShading
    imagesc(z,Y,P);
    set(gca,'Layer','top'); % show grid lines on top
end;

% Plot mean and std
plot(z, m, 'color',col_d,'LineWidth',lineWidthMean);
plot(z, m+n_std*std_V, z, m-n_std*std_V,'color',col_d,'LineWidth',lineWidthStd);
hold on;
if ~flag_noSamples
    for iS=1:N_sample
        plot(z,s_prior(:,iS),'--','color',col_d,'LineWidth',lineWidthStd);
    end;
end;

% Axis settings
extraLim = 0.03;
set(gca,'xlim',[z(1)-extraLim, z(end)+extraLim]);
set(gca,'ylim',[min_std_V, max_std_V]);
%set(gca,'ylim',[ min([min(s_prior), -multSTD*max_std_V])-extraLim, ...
%    max([max(s_prior), multSTD*max_std_V])+extraLim ]);

drawnow;
