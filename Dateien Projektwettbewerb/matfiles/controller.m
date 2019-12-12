function [U] = controller(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function [U] = controller(X)
%
% controller for the single-track model
%
% inputs: x (x position), y (y position), v (velocity), beta
% (side slip angle), psi (yaw angle), omega (yaw rate), x_dot (longitudinal
% velocity), y_dot (lateral velocity), psi_dot (yaw rate (redundant)), 
% varphi_dot (wheel rotary frequency)
%
% external inputs (from 'racetrack.mat'): t_r_x (x coordinate of right 
% racetrack boundary), t_r_y (y coordinate of right racetrack boundary),
% t_l_x (x coordinate of left racetrack boundary), t_l_y (y coordinate of
% left racetrack boundary)
%
% outputs: delta (steering angle ), G (gear 1 ... 5), F_b (braking
% force), zeta (braking force distribution), phi (gas pedal position)
%
% files requested: racetrack.mat
%
% This file is for use within the "Project Competition" of the "Concepts of
% Automatic Control" course at the University of Stuttgart, held by F.
% Allgoewer.
%
% prepared by J. M. Montenbruck, Dec. 2013 
% mailto:jan-maximilian.montenbruck@ist.uni-stuttgart.de
%
% written by *STUDENT*, *DATE*
% mailto:*MAILADDRESS*


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% state vector
x=X(1); % x position
y=X(2); % y position
v=X(3); % velocity (strictly positive)
beta=X(4); % side slip angle
psi=X(5); % yaw angle
omega=X(6); % yaw rate
x_dot=X(7); % longitudinal velocity
y_dot=X(8); % lateral velocity
psi_dot=X(9); % yaw rate (redundant)
varphi_dot=X(10); % wheel rotary frequency (strictly positive)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% racetrack
load('racetrack.mat','t_r'); % load right  boundary from *.mat file
load('racetrack.mat','t_l'); % load left boundary from *.mat file
t_r_x=t_r(:,1); % x coordinate of right racetrack boundary
t_r_y=t_r(:,2); % y coordinate of right racetrack boundary
t_l_x=t_l(:,1); % x coordinate of left racetrack boundary
t_l_y=t_l(:,2); % y coordinate of left racetrack boundary
dt = 0.1;
%% 
d = @(z)deal(0,0);
sigmaw = 0;
% estModel = MotionModelGP_SingleTrack(d,sigmaw);
estModel = MotionModelGP_TrueSingleTrack(d,sigmaw);

n = estModel.n;
m = estModel.m;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NONLINEAR MPC CONTROLLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define cost function
N = 10;     % prediction horizon
Q = [100,0 ;0,100];
Qf= [100,0;0,100];
R = eye(5);
Ck = [eye(2), zeros(2,8)];
fo   = @(t,mu_x,var_x,u,r) (Ck*mu_x-r(t))'*Q *(Ck*mu_x-r(t)) + u'*R*u;  % cost function
fend = @(t,mu_x,var_x,r)   (Ck*mu_x-r(t))'*Qf*(Ck*mu_x-r(t));          % end cost function
f    = @(mu_xk,var_xk,u) estModel.xkp1(mu_xk, var_xk, u, dt);
h    = @(x,u) []; % @(x,u) 0;  % h(x)==0
g    = @(x,u) []; % @(x,u) 0;  % g(x)<=0

% mpc = NMPC(fo, fend, f, d_GP, Bd, Bz, N, sigmaw, h, g, n, m, ne, dt);
mpc = NMPC(f, h, g, n, m, fo, fend, N, dt);
mpc.tol     = 1e-3;
mpc.maxiter = 30;

% NMPC controller
% 

% calculate trajectory center line
t_c = (t_r + t_l)/2;
% add a few points at end
t_c = [t_c; t_c(2:20,:)];
% find closest trajectory point w.r.t. the vehicle
[~,idx] = min( pdist2(X(1:2)',t_c,'seuclidean',[1 1].^0.5).^2 );
% set target as 3 poins ahead
idx_target = idx + 8;
% loop around when track is over
idx_target = mod(idx_target, size(t_c,1));

r = @(t) [t_c(idx_target,:)]';
U = mpc.optimize(X, 0, r);

X(1:3)
r(3234)


%% stop condition
if X(1)<0 & X(1)>-5 & X(2)<0 & X(2)>-1
    return
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% U=[delta G F_b zeta phi]; % input vector




end

