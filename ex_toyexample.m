clear vars; close all; clc;

dt = 0.01;

%% Dynamic Model
%------------------------------------------------------------------
%   dot_x(t) = f_true(x(t),u(t)) + Bd*(d_true(t) + w(t)),    wk=N(0,sigmaw^2)
%   
%   x = [x1]
%   u = [u1]
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


%% Controller

% nominal model
f = @(x,u) A*x + B*u;
% discretize nominal model
fd = @(x,u,dt) x + dt*f(x,u);

% LQR controller
Ak = eye(n)+dt*A;
Bk = B*dt;
K = place(Ak,Bk,0.9);

% Prefilter
Kr = inv((eye(n)-Ak+Bk*K)\Bk);

eig(Ak-Bk*K)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO:
%       - Implement MPC controller - See MPC.m class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Gaussian Process

% GP hyperparameters
sigmaf  = 0.1;              % output variance (std)
lambda  = diag([2,10].^2);   % length scale
sigman  = sigmaw*sqrt(dt);  % stddev of measurement noise
maxsize = 100;              % maximum number of points in the dictionary

% create GP object
gp = GP(sigmaf, sigman, lambda, maxsize);


%% Simulate

% define input
% r = @(t) -3;
% r = @(t) 1*sin(10*t);
r = @(t) 5*sin(5*t) + 5*sin(15*t) + 10*exp(-t);
nr = 1; % dimension of r(t)

% simulation time
tf = 10;

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
    out.u(:,i) = Kr*out.r(:,i) - K*out.xhat(:,i);
    
    % simulate real model
    out.x(:,i+1) = fd_true(out.x(:,i),out.u(:,i),dt,true);
    
    % measure data
    out.xhat(:,i+1) = out.x(:,i+1); % perfect observer
    
    % add data to GP model
    out.xnom(:,i+1) = fd(out.xhat(:,i),out.u(:,i),dt);
    if mod(i-1,7)==0
        ddata = Bd \ (out.xhat(:,i+1) - out.xnom(:,i+1));
        gp.add( [out.xnom(:,i);out.u(:,i)], ddata );
    end
    
%     ddata
%     gptrue( [out.xhat(:,i),out.u(:,i)] )
%     d_true( out.xhat(:,i),out.u(:,i) )*dt
    
end


%% Evaluate results
close all;

% true GP function that is ment to be learned
gptrue = @(x) Bd \ ( fd_true(x(1),x(2),dt,false) - fd(x(1),x(2),dt) );

% plot prediction bias and variance
gp.plot2d( gptrue )
       
% plot reference and state signal
figure; hold on; grid on;
plot(out.t(1:end-1), out.r, 'DisplayName', 'r(t)')
plot(out.t, out.x, 'DisplayName', 'x(t)')
legend;




