function simulate_controlled_singletrack(t_f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function simulate_controlled_singletrack(t_f)
%
% integrates the controlled single-track model until time t_f
%
% input: t_f (simulation time)
%
% files requested: racetrack.m ; singletrack.m ; ode1.m ; plot_racetrack.m
%
% plots built: racetrack
%
% This file is for use within the "Project Competition" of the "Concepts of
% Automatic Control" course at the University of Stuttgart, held by F.
% Allgoewer.
%
% written by J. M. Montenbruck, Dec. 2013 
% mailto:jan-maximilian.montenbruck@ist.uni-stuttgart.de

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
racetrack % builds the racetrack and saves it as racetrack.mat


dt = 0.05;


% -------------------------------------------------------------------------
% NONLINEAR MPC CONTROLLER
% define cost function
N = 20;     % prediction horizon
Q = 100;
Qf= 100;
R = 1;
f    = @(t,x,u) singletracknominal(x,u,dt);     % nominal model
fo   = @(t,x,u,e,r) (x-r(t))'*Q *(x-r(t)) + R*u^2;  % cost function
fend = @(t,x,e,r)   (x-r(t))'*Qf*(x-r(t));  % end cost function
h    = []; % @(t,x,u,e) 0;  % h(x)==0
g    = []; % @(t,x,u,e) 0;  % g(x)<=0
ne   = 0;

mpc = NMPC(fo, fend, f, d_gp, Bd, N, sigmaw, h, g, n, m, ne, dt);
mpc.tol     = 1e-3;
mpc.maxiter = 30;
% -------------------------------------------------------------------------




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTEGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_0=[-2.5;0;0;0;pi/2;0;0;0;0;0]; % initial value for integration
Y=ode1(@singletrack,0:dt:t_f,X_0); % integrate with step zise 0.001

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EVALUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_racetrack % plots the racetrack and your result
end