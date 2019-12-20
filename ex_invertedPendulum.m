%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -
%------------------------------------------------------------------

%------------------------------------------------------------------
%   1D Toy example:
%
%       - Simulate the GP learning of the nonlinear part of the plant
%       dynamics
%       - System is being currently controlled with a state feedback control
%------------------------------------------------------------------

clear all; close all; clc;

dt = 0.1;  % simulation timestep size
tf = 7;     % simulation time

% inverted pendulum parameters
Mc = 5;
Mp = 2;
b = 0.1;
I = 0.6;
l = 3;



%% True Dynamics Model
%--------------------------------------------------------------------------
%   xk+1 = fd(xk,uk) + Bd * ( d(zk) + w ),    
%
%       where: zk = Bz*xk,
%              d ~ N(mean_d(zk),var_d(zk))
%              w ~ N(0,var_w)
%------------------------------------------------------------------

% define model (mean and variance) for true disturbance
% mu_d  = @(z) 1 * mvnpdf(z',[0,0], eye(2)*0.1);
mu_d  = @(z) 0.1 * z(1) - 0.01*z(2) + deg2rad(3);
var_d = @(z) 0;
d_true  = @(z) deal(mu_d(z),var_d(z));
% true measurement noise
var_w = 1e-8;
% create true dynamics model
trueModel = MotionModelGP_InvertedPendulum(Mc, Mp, b, I, l, d_true, var_w);


%% Create Estimation Model and Nominal Model

% define model (mean and variance) for estimated disturbance
% GP hyperparameters
var_f   = 0.01;                     % output variance (std)
M       = diag([1e-1,1e-1].^2);     % length scale
var_n   = var_w;                    % measurement noise variance
maxsize = 100;                      % maximum number of points in the dictionary
% create GP object
d_GP = GP(var_f, var_n, M, maxsize);


% create estimation dynamics model (disturbance is the Gaussian Process GP)
estModel = MotionModelGP_InvertedPendulum(Mc, Mp, b, I, l, @d_GP.eval, var_w);

% create nominal dynamics model (no disturbance)
nomModel = MotionModelGP_InvertedPendulum(Mc, Mp, b, I, l, @(z)deal(0,0), 0); 


%% Controller

n = estModel.n;
m = estModel.m;
ne = 0;

% -------------------------------------------------------------------------
% LQR CONTROLLER
[A,B] = estModel.linearize();
Ak = eye(n)+dt*A;
Bk = B*dt;
Ck=[0 1 0 0; 0 0 1 0; 0 0 0 1];
Q = 1e3*eye(4);
R = 1;
[~,~,K] = dare(Ak,Bk,Q,R);
% Prefilter
Kr = pinv(Ck/(eye(n)-Ak+Bk*K)*Bk);
% check eigenvalues
eig(Ak-Bk*K);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% NONLINEAR MPC CONTROLLER
% define cost function
N = 10;     % prediction horizon
Q = diag([1e-1 1e5 1]);
Qf= diag([1e-1 1e5 1]);
R = 1;
Ck = [0 1 0 0; 0 0 1 0; 0 0 0 1];
fo   = @(t,mu_x,var_x,u,e,r) (Ck*mu_x-r(t))'*Q *(Ck*mu_x-r(t)) + R*u^2;  % cost function
fend = @(t,mu_x,var_x,e,r)   (Ck*mu_x-r(t))'*Qf*(Ck*mu_x-r(t));          % end cost function
f    = @(mu_xk,var_xk,u) estModel.xkp1(mu_xk, var_xk, u, dt);
h    = @(x,u,e) []; % @(x,u) 0;  % h(x)==0
g    = @(x,u,e) []; % @(x,u) 0;  % g(x)<=0

u_lb = [];
u_ub = [];

mpc = NMPC (f, h, g, u_lb, u_ub, n, m, ne, fo, fend, N, dt); % (f, h, g, n, m, fo, fend, N, dt);
mpc.tol     = 1e-3;
mpc.maxiter = 30;
% -------------------------------------------------------------------------

% TEST NMPC
x0 = [0 0 0.1 0]';
t0 = 0;
r  = @(t) [0 0 0]';    % desired trajectory
[x0_opt, u_opt, e_opt] = mpc.optimize(x0, t0, r );



%% Simulate

% define input
r = @(t) [0 0 0]';
% r = @(t) 1*sin(10*t);
% r = @(t) 2*sin(5*t) + 2*sin(15*t) + 6*exp(-t) - 4 ;
% r = @(t) 4*sin(5*t) + 4*sin(15*t);
nr = size(r(0),1); % dimension of r(t)

% initial state
x0 = [0,0,deg2rad(5),0]';

% initialize variables to store simulation results
out.t = 0:dt:tf;
out.x = [x0 zeros(n,length(out.t)-1)];
out.xhat = [x0 zeros(n,length(out.t)-1)];
out.xnom = [x0 zeros(n,length(out.t)-1)];
out.u = zeros(m,length(out.t)-1);
out.r = zeros(nr,length(out.t)-1);


d_GP.isActive = false;

for k = 1:numel(out.t)-1
    disp(out.t(k))
    
    % ---------------------------------------------------------------------
    % Read new reference
    % ---------------------------------------------------------------------
    out.r(:,k) = r(out.t(k));
    
    % ---------------------------------------------------------------------
    % LQR controller
    % ---------------------------------------------------------------------
    % % out.u(:,i) = Kr*out.r(:,i) - K*out.xhat(:,i);
    
    % ---------------------------------------------------------------------
    % NPMC controller
    % ---------------------------------------------------------------------
    [x0_opt, u_opt, e_opt] = mpc.optimize(out.xhat(:,k), out.t(k), r);
    out.u(:,k) = u_opt(:,1);
    
    
    % ---------------------------------------------------------------------
    % simulate real model
    % ---------------------------------------------------------------------
    [mu_xkp1,var_xkp1] = trueModel.xkp1(out.x(:,k),zeros(trueModel.n),out.u(:,k),dt);
    out.x(:,k+1) = mvnrnd(mu_xkp1, var_xkp1, 1)';
    
    
    % ---------------------------------------------------------------------
    % measure data
    % ---------------------------------------------------------------------
    out.xhat(:,k+1) = out.x(:,k+1); % perfect observer
    
    
    % ---------------------------------------------------------------------
    % calculate nominal model
    % ---------------------------------------------------------------------
    out.xnom(:,k+1) = nomModel.xkp1(out.xhat(:,k),zeros(nomModel.n),out.u(:,k),dt);
    
    
    % ---------------------------------------------------------------------
    % add data to GP model
    % ---------------------------------------------------------------------
    if mod(k-1,1)==0
        % calculate disturbance (error between measured and nominal)
        d_est = estModel.Bd \ (out.xhat(:,k+1) - out.xnom(:,k+1));
        % select subset of coordinates that will be used in GP prediction
        zhat = estModel.Bz * out.xhat(:,k);
        % add data point to the GP dictionary
        d_GP.add(zhat,d_est);
    end
    
    if d_GP.N > 20 && out.t(k) > 3
        d_GP.updateModel();
        d_GP.isActive = true;
    end
    
    % check if these values are the same:
    % d_est == mu_d(zhat) == ([mud,~]=trueModel.d(zhat); mud*dt)
    
end


%% Evaluate results
close all;

d_GP.isActive = true;

% plot reference and state signal
figure('Position',[-1836 535 560 420]); 
subplot(2,1,1); hold on; grid on;
plot(out.t(1:end-1), out.r, 'DisplayName', 'r(t)')
plot(out.t, out.x(3,:), 'DisplayName', 'x(t) [rad]')
ylabel('[rad]');
legend;
subplot(2,1,2); hold on; grid on;
plot(out.t(1:end-1), out.u, 'DisplayName', 'u(t)')
legend;

% true GP function that is ment to be learned
Bz = trueModel.Bz;
Bd = trueModel.Bd;

% define the true expected disturbance model
% z = [0;0.1];
% gptrue = @(z) mu_d(z);
gptrue = @(z) trueModel.Bd'*( trueModel.xkp1(trueModel.Bz'*z,zeros(trueModel.n),0,dt)...
                             -nomModel.xkp1(trueModel.Bz'*z,zeros(nomModel.n),0,dt)  );

% plot prediction bias and variance
d_GP.plot2d( gptrue )
       





