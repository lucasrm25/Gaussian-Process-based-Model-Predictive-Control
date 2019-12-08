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
%              w ~ N(0,sigmaw)
%------------------------------------------------------------------

% define model (mean and variance) for true disturbance
% mu_d  = @(z) 1 * mvnpdf(z',[0,0], eye(2)*0.1);
mu_d  = @(z) 0.1 * z(1) - 0.01*z(2) + deg2rad(3);
var_d = @(z) 0;
d_true  = @(z) deal(mu_d(z),var_d(z));
% true measurement noise
sigmaw = 1e-8;
% create true dynamics model
trueModel = invertedPendulumModel(Mc, Mp, b, I, l, d_true, sigmaw);



%% Create Estimation Model and Nominal Model

% define model (mean and variance) for estimated disturbance
% GP hyperparameters
sigmaf2 = 0.01;         % output variance (std)
M       = diag([1e-1,1e-1].^2);   % length scale
sigman2 = 1e-5;         % measurement noise variance
maxsize = 100;          % maximum number of points in the dictionary
% create GP object
d_GP = GP(sigmaf2, sigman2, M, maxsize);


% create estimation dynamics model (disturbance is the Gaussian Process GP)
estModel = invertedPendulumModel(Mc, Mp, b, I, l, @d_GP.eval, sigmaw);

% create nominal dynamics model (no disturbance)
nomModel = invertedPendulumModel(Mc, Mp, b, I, l, @(z)deal(0,0), 0); 


%% Controller

n = estModel.n;
m = estModel.m;

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
fo   = @(t,x,u,r) (Ck*x-r(t))'*Q *(Ck*x-r(t)) + R*u^2;  % cost function
fend = @(t,x,r)   (Ck*x-r(t))'*Qf*(Ck*x-r(t));            % end cost function
f    = @(x,u) estModel.xkp1(x,u,dt);
h    = @(x,u) []; % @(x,u) 0;  % h(x)==0
g    = @(x,u) []; % @(x,u) 0;  % g(x)<=0

% mpc = NMPC(fo, fend, f, d_GP, Bd, Bz, N, sigmaw, h, g, n, m, ne, dt);
mpc = NMPC(f, h, g, n, m, fo, fend, N, dt);
mpc.tol     = 1e-3;
mpc.maxiter = 30;
% -------------------------------------------------------------------------

% TEST NMPC
% x0 = 10;
% t0 = 0;
% r  = @(t)2;    % desired trajectory
% u0 = mpc.optimize(x0, t0, r );



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

for i = 1:numel(out.t)-1
    disp(out.t(i))
    
    % read new reference
    out.r(:,i) = r(out.t(i));
    
    % calculate control input
    % LQR controller
    % % out.u(:,i) = Kr*out.r(:,i) - K*out.xhat(:,i);
    % NMPC controller
    out.u(:,i) = mpc.optimize(out.xhat(:,i), out.t(i), r);
    
    % simulate real model
    [mu_xkp1,var_xkp1] = trueModel.xkp1(out.x(:,i),out.u(:,i),dt);
    out.x(:,i+1) = mvnrnd(mu_xkp1, var_xkp1, 1)';
    
    % measure data
    out.xhat(:,i+1) = out.x(:,i+1); % perfect observer
    
    % calculate nominal model
    out.xnom(:,i+1) = nomModel.xkp1(out.xhat(:,i),out.u(:,i),dt);
    
    % add data to GP model
    if mod(i-1,1)==0
        % calculate disturbance (error between measured and nominal)
        d_est = estModel.Bd \ (out.xhat(:,i+1) - out.xnom(:,i+1));
        % select subset of coordinates that will be used in GP prediction
        zhat = estModel.Bz * out.xhat(:,i);
        % add data point to the GP dictionary
        d_GP.add(zhat,d_est);
    end
    
    if d_GP.N > 20 && out.t(i) > 3
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
gptrue = @(z) trueModel.Bd'*(trueModel.xkp1(trueModel.Bz'*z,0,dt) - nomModel.xkp1(trueModel.Bz'*z,0,dt));

% plot prediction bias and variance
d_GP.plot2d( gptrue )
       





