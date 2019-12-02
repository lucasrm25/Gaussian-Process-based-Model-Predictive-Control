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

clear vars; close all; clc;

dt = 0.05;  % simulation timestep size
tf = 2;     % simulation time

%% Dynamic Model
%------------------------------------------------------------------
%   dot_x(t) = f_true(x(t),u(t)) + Bd*(d_true(t) + w(t)),    wk=N(0,sigmaw^2)
%   
%   x = [x1]        1D state
%   u = [u1]        1D input
%------------------------------------------------------------------

% true model
A = -0.1;
B = 3;
Bd = 1;

% disturbance noise stddev - in continuous time
sigmaw = 0.01/sqrt(dt);

n = size(A,2);
m = size(B,2);
md = size(Bd,2);

% true disturbance model
d_true = @(x,u) 1/dt*[20 -20]*[mvnpdf([x,u],[2,2],diag([2,20])) mvnpdf([x,u],[-2,-2],diag([2,20]))]';

% true model
f_true = @(x,u) A*x + B*u + Bd*(d_true(x,u));

% dicretize true model - ODE1 Euler integration
fd_true = @(x,u,dt,inclnoise) x + dt*f_true(x,u) + inclnoise*sqrt(dt)*Bd*sigmaw*randn(md);



%% Gaussian Process

% GP hyperparameters
sigmaf  = 0.1;              % output variance (std)
lambda  = diag([2,10].^2);  % length scale
sigman  = sigmaw*sqrt(dt);  % stddev of measurement noise
maxsize = 100;              % maximum number of points in the dictionary

% create GP object
d_gp = GP(sigmaf, sigman, lambda, maxsize);



%% Controller

% -------------------------------------------------------------------------
% DEFINE NOMINAL MODEL
% nominal continuous time model
f = @(x,u) A*x + B*u;
% discretize nominal model
fd = @(x,u,dt) x + dt*f(x,u);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% LQR CONTROLLER
Ak = eye(n)+dt*A;
Bk = B*dt;
K = place(Ak,Bk,0.9);
% Prefilter
Kr = inv((eye(n)-Ak+Bk*K)\Bk);
% check eigenvalues
eig(Ak-Bk*K);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% NONLINEAR MPC CONTROLLER
% define cost function
N = 20;     % prediction horizon
Q = 100;
Qf= 100;
R = 1;
f    = @(t,x,u) fd(x,u,dt);
fo   = @(t,x,u,e,r) (x-r(t))'*Q *(x-r(t)) + R*u^2;
fend = @(t,x,e,r)   (x-r(t))'*Qf*(x-r(t));
h    = []; % @(t,x,u,e) 0;
g    = []; % @(t,x,u,e) 0;
ne   = 0;

mpc = NMPC(fo, fend, f, d_gp, Bd, N, sigmaw, h, g, n, m, ne, dt);
mpc.tol     = 1e-3;
mpc.maxiter = 30;

x0 = 10;
e0 = 0;
t0 = 0;
r  = @(t)2;    % desired trajectory
u0 = mpc.optimize(x0, e0, t0, r );
% -------------------------------------------------------------------------


%% Simulate

% define input
% r = @(t) -3;
% r = @(t) 1*sin(10*t);
r = @(t) 2*sin(5*t) + 2*sin(15*t) + 6*exp(-t) - 4 ;
nr = 1; % dimension of r(t)

% initial state
x0 = 0;

% initialize variables to store simulation results
out.t = 0:dt:tf;
out.x = [x0 zeros(n,length(out.t)-1)];
out.xhat = [x0 zeros(n,length(out.t)-1)];
out.xnom = [x0 zeros(n,length(out.t)-1)];
out.u = zeros(m,length(out.t)-1);
out.r = zeros(nr,length(out.t)-1);


for i = 1:numel(out.t)-1
    disp(out.t(i))
    
    % read new reference
    out.r(:,i) = r(out.t(i));
    
    % calculate control input
    % out.u(:,i) = Kr*out.r(:,i) - K*out.xhat(:,i);
    x0 = out.xhat(:,i);
    e0 = 0;
    t0 = out.t(i);
    out.u(:,i) = mpc.optimize(x0, e0, t0, r);
    
    % simulate real model
    out.x(:,i+1) = fd_true(out.x(:,i),out.u(:,i),dt,true);
    
    % measure data
    out.xhat(:,i+1) = out.x(:,i+1); % perfect observer
    
    
    % add data to GP model
    out.xnom(:,i+1) = fd(out.xhat(:,i),out.u(:,i),dt);
    if mod(i-1,2)==0
        % calculate disturbance (error between measured and nominal)
        disturb = Bd \ (out.xhat(:,i+1) - out.xnom(:,i+1));
        % select subset of coordinates that will be used in GP prediction
        zhat = Bd*out.xhat(:,i);
        % add data point to the GP dictionary
        mpc.d_gp.add( [zhat;out.u(:,i)], disturb );
    end
    
    if mpc.d_gp.N > 50
        mpc.activateGP();
    end
    
    % check if these tree values are the same:
    %     ddata
    %     gptrue( [out.xhat(:,i),out.u(:,i)] )
    %     d_true( out.xhat(:,i),out.u(:,i) )*dt
    
end


%% Evaluate results
close all;

% true GP function that is ment to be learned
gptrue = @(x) Bd \ ( fd_true(x(1),x(2),dt,false) - fd(x(1),x(2),dt) );

% plot prediction bias and variance
mpc.d_gp.plot2d( gptrue )
       
% plot reference and state signal
figure; 
subplot(2,1,1); hold on; grid on;
plot(out.t(1:end-1), out.r, 'DisplayName', 'r(t)')
plot(out.t, out.x, 'DisplayName', 'x(t)')
legend;
subplot(2,1,2); hold on; grid on;
plot(out.t(1:end-1), out.u, 'DisplayName', 'u(t)')
legend;




