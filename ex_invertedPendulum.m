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

dt = 0.05;  % simulation timestep size
tf = 2.5;     % simulation time

% inverted pendulum parameters
Mc = 5;
Mp = 2;
b = 0.1;
I = 0.6;
l = 3;



%% True Dynamics Model
%------------------------------------------------------------------
%   dot_x(t) = f_true(x(t),u(t)) + Bd*(d(z(t))),    d~N(mu_d(z),var_d(z))
%------------------------------------------------------------------

% define model (mean and variance) for true disturbance
% mu_d  = @(z) 1 * mvnpdf(z',[0,0], eye(2)*0.1);
mu_d  = @(z) 0.5 * z(1)^0.5 - 0.1*z(2);
var_d = @(z) 0*1e-5;
d_true  = @(z) deal(mu_d(z),var_d(z));

% create true dynamics model
trueModel = invertedPendulum(Mc, Mp, b, I, l, d_true);



%% Create Estimation Model and Nominal Model

% define model (mean and variance) for estimated disturbance
% GP hyperparameters
sigmaf2 = 0.01;         % output variance (std)
lambda  = diag([1e-1,1e-1].^2);   % length scale
sigman2 = 1e-5;         % measurement noise variance
maxsize = 100;          % maximum number of points in the dictionary
% create GP object
d_GP = GP(sigmaf2, sigman2, lambda, maxsize);


% create estimation dynamics model (disturbance is the Gaussian Process GP)
estModel = invertedPendulum(Mc, Mp, b, I, l, @d_GP.eval);

% create nominal dynamics model (no disturbance)
nomModel = invertedPendulum(Mc, Mp, b, I, l, @(z) deal(0,0) ); 


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
Q = diag([1e-1 1e4 1]);
Qf= diag([1e-1 1e4 1]);
R = 1;
Ck = [0 1 0 0; 0 0 1 0; 0 0 0 1];
fo   = @(t,x,u,r) (Ck*x-r(t))'*Q *(Ck*x-r(t)) + R*u^2;  % cost function
fend = @(t,x,r)   (Ck*x-r(t))'*Qf*(Ck*x-r(t));            % end cost function
f    = @(x,u) estModel.fd(x,u,dt);
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
    x0 = out.xhat(:,i);
    t0 = out.t(i);
    out.u(:,i) = mpc.optimize(x0, t0, r);
    
    % simulate real model
    [mu_xkp1,var_xkp1] = trueModel.fd(out.x(:,i),out.u(:,i),dt);
    out.x(:,i+1) = mvnrnd(mu_xkp1, var_xkp1, 1)';
    
    % measure data
    out.xhat(:,i+1) = out.x(:,i+1); % perfect observer
    
    % calculate nominal model
    out.xnom(:,i+1) = nomModel.fd(out.xhat(:,i),out.u(:,i),dt);
    
    % add data to GP model
    if mod(i-1,1)==0
        % calculate disturbance (error between measured and nominal)
        d_est = estModel.Bd \ (out.xhat(:,i+1) - out.xnom(:,i+1));
        % select subset of coordinates that will be used in GP prediction
        zhat = estModel.Bz * out.xhat(:,i);
        % add data point to the GP dictionary
        d_GP.add(zhat,d_est);
    end
    
    if d_GP.N > 20 && out.t(i) > 2
        d_GP.isActive = true;
    end
    
    % check if these values are the same:
    % d_est == mu_d(zhat)*dt == ([mud,~]=trueModel.d(zhat); mud*dt)
end


%% Evaluate results
close all;

d_GP.isActive = true;

% plot reference and state signal
figure; 
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
% gptrue = @(z) mu_d(z)*dt; % also works
gptrue = @(z) trueModel.Bd'*(trueModel.fd(trueModel.Bz'*z,0,dt) - nomModel.fd(trueModel.Bz'*z,0,dt));

% plot prediction bias and variance
d_GP.plot2d( gptrue )
       





