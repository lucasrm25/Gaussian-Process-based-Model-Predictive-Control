%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -

%   1D Toy example:
%
%       - Simulate the GP learning of the nonlinear part of the plant
%       dynamics
%       - System is being currently controlled with a state feedback control
%------------------------------------------------------------------

clear all; close all; clc;



%--------------------------------------------------------------------------
%   Quick Access Simulation and controller parameters
%------------------------------------------------------------------
dt = 0.1;       % simulation timestep size
tf = 7;         % simulation time
maxiter = 15;   % max NMPC iterations per time step
N = 10;         % NMPC prediction horizon


useParallel = false;


lookahead = dt*N;
fprintf('\nPrediction lookahead: %.1f [s]\n',lookahead);


% inverted pendulum parameters
Mc = 5;
Mp = 2;
b = 0.1;
I = 0.6;
l = 3;
g = 9.81;


%% True Dynamics Model
%--------------------------------------------------------------------------
%   xk+1 = fd_true(xk,uk) + Bd * ( w ),    
%
%       where: w ~ N(0,var_w)
%------------------------------------------------------------------

% define noise for true disturbance
var_w = 1e-8;

% create true dynamics model
trueModel = MotionModelGP_InvPendulum_deffect(Mc, Mp, b, I, l, [], var_w);

%% Create Estimation Model and Nominal Model

% -------------------------------------------------------------------------
%  Create nominal model (no disturbance):  
%       xk+1 = fd_nom(xk,uk)
% -------------------------------------------------------------------------

% create nominal dynamics model (no disturbance)
nomModel = MotionModelGP_InvPendulum_nominal(Mc, Mp, b, I, l, [], []); 


% -------------------------------------------------------------------------
%  Create adaptive dynamics model 
%  (unmodeled dynamics will be estimated by Gaussian Process GP)
%       xk+1 = fd_nom(xk,uk) + Bd * ( d_GP(zk) + w )
% -------------------------------------------------------------------------

% GP input dimension
gp_n = MotionModelGP_InvPendulum_nominal.nz;
% GP output dimension
gp_p = MotionModelGP_InvPendulum_nominal.nd;

% GP hyperparameters
var_f   = 0.01;                     % output variance
M       = diag([1e-1,1e-1].^2);     % length scale
var_n   = var_w;                    % measurement noise variance
maxsize = 100;                      % maximum number of points in the dictionary

% create GP object
d_GP = GP(gp_n, gp_p, var_f, var_n, M, maxsize);

% create estimation dynamics model (disturbance is the Gaussian Process GP)
estModel = MotionModelGP_InvPendulum_nominal(Mc, Mp, b, I, l, @d_GP.eval, var_w);



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
Q = diag([1e-1 1e5 1e0]);
Qf= diag([1e-1 1e5 1e0]);
R = 10;
Ck = [0 1 0 0; 0 0 1 0; 0 0 0 1];
fo   = @(t,mu_x,var_x,u,e,r) (Ck*mu_x-r(t))'*Q *(Ck*mu_x-r(t)) + R*u^2;  % cost function
fend = @(t,mu_x,var_x,e,r)   (Ck*mu_x-r(t))'*Qf*(Ck*mu_x-r(t));          % end cost function
f    = @(mu_xk,var_xk,u) estModel.xkp1(mu_xk, var_xk, u, dt);
h    = @(x,u,e) []; % @(x,u) 0;  % h(x)==0
g    = @(x,u,e) []; % @(x,u) 0;  % g(x)<=0

u_lb = [];
u_ub = [];

mpc = NMPC (f, h, g, u_lb, u_ub, n, m, ne, fo, fend, N, dt);
mpc.tol     = 1e-3;
mpc.maxiter = maxiter;
% -------------------------------------------------------------------------


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
out.t    = 0:dt:tf;
out.x    = [x0 nan(n,length(out.t)-1)];
out.xhat = [x0 nan(n,length(out.t)-1)];
out.xnom = [x0 nan(n,length(out.t)-1)];
out.u    = nan(m,length(out.t)-1);
out.r    = nan(nr,length(out.t)-1);


d_GP.isActive = false;


ki = 1;
% ki = 40;
% mpc.uguess = out.u(:,ki);

for k = ki:numel(out.t)-1
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
    [u_opt, e_opt] = mpc.optimize(out.xhat(:,k), out.t(k), r, useParallel);
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
    % Safety
    % ---------------------------------------------------------------------
    if abs(out.xhat(3,k+1)) > deg2rad(60)
        error('Pole is completely unstable. theta = %.f[deg]... aborting',rad2deg(out.xhat(3,k+1)));
    end
    
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
        zhat = [ estModel.Bz_x * out.xhat(:,k); estModel.Bz_u * out.u(:,k) ];
        % add data point to the GP dictionary
        d_GP.add(zhat,d_est);
    end
    
    if d_GP.N > 20 && out.t(k) > 3
        d_GP.updateModel();
        d_GP.isActive = true;
    end
    
    % check if these values are the same:
    % d_est == mu_d(zhat) == [mud,~]=trueModel.d(zhat)
    
end




return




%% Optimize GP hyperparameters ??? (Offline procedure, after simulation)

d_GP.setHyperParameters( M, var_f, var_n )
% d_GP.optimizeHyperParams('ga');
d_GP.optimizeHyperParams('fmincon');

d_GP.M
d_GP.var_f
d_GP.var_n


%% Evaluate results
close all;

% plot reference and state signal
figure('Color','w','Position',[-1836 535 560 420]); 
subplot(2,1,1); hold on; grid on;
% plot(out.t(1:end-1), out.r, 'DisplayName', 'r(t)')
plot(out.t, out.x(3,:), 'DisplayName', 'x(t) [rad]')
ylabel('Pole angle \theta [rad]');
xlabel('time [s]')

subplot(2,1,2); hold on; grid on;
plot(out.t(1:end-1), out.u, 'DisplayName', 'u(t)')
ylabel('Force on the carriage F [N]');
xlabel('time [s]')
% legend;

% true GP function that is meant to be learned
Bz_x = trueModel.Bz_x;
Bz_u = trueModel.Bz_u;
Bd = trueModel.Bd;
n = trueModel.n;

% define the true expected disturbance model
% z = [0;0.1];
gptrue = @(z) Bd'*( trueModel.xkp1(Bz_x'*z, zeros(n), 0, dt)...
                   - nomModel.xkp1(Bz_x'*z, zeros(n), 0, dt)  );

% plot prediction bias and variance
d_GP.plot2d( gptrue )

%% animation of inverse pendulum

% animation of inverse pendulum
drawpendulum(out.t,out.x,Mc,Mp,g,l)     


%% Analyse learning
% ---------------------------------------------------------------------
% Check how the GP reduces the prediction error
% ---------------------------------------------------------------------

% d_GP.optimizeHyperParams('fmincon')
% d_GP.optimizeHyperParams('ga')


k = find(~isnan(out.xhat(1,:)), 1, 'last' ) - 1;

% prediction error without GP
% predErrorNOgp = estModel.Bd\(out.xhat - out.xnom);
predErrorNOgp = estModel.Bd\(out.xhat(:,1:k-1) - out.xnom(:,1:k-1));


% prediction error with trained GP
zhat = estModel.z( out.xhat(:,1:k-1), out.u(:,1:k-1) )
dgp = d_GP.eval(zhat,true);
predErrorWITHgp = estModel.Bd\( out.xhat(:,2:k) - (out.xnom(:,2:k) + estModel.Bd*dgp) );


disp('Prediction mean squared error without GP:')
disp( mean(predErrorNOgp(:,all(~isnan(predErrorNOgp))).^2 ,2) )
disp('Prediction mean squared error with trained GP:')
disp( mean(predErrorWITHgp(:,all(~isnan(predErrorWITHgp))).^2 ,2) )



% Visualize error
figure('Color','w'); hold on; grid on;
subplot(1,2,1)
plot( predErrorNOgp' )
subplot(1,2,2)
hist(predErrorNOgp')
sgtitle('Prediction error - without GP')


figure('Color','w'); hold on; grid on;
subplot(1,2,1)
plot( predErrorWITHgp' )
subplot(1,2,2)
hist(predErrorWITHgp')
sgtitle('Prediction error - with GP')



