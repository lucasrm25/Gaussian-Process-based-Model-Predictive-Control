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
sigmaw = 0.1/dt;

n = size(A,2);
m = size(B,2);
md = size(Bd,2);

% true disturbance model
d_true = @(x,u) -3*[0.01 -1 1 0.5]*[(x+u)^3 exp(-10*(u+0.5)^2) exp(-10*(x-0.5)^2) sin(2*x)]';
% true model
f_true = @(x,u,inclnoise) A*x + B*u + Bd*(d_true(x,u) + inclnoise*sigmaw*randn(md));

% dicretize true model - ODE1 Euler integration
fd_true = @(x,u,dt,inclnoise) x + dt*f_true(x,u,inclnoise);


%% Controller

% nominal model
f = @(x,u) A*x + B*u;
% discretize nominal model
fd = @(x,u,dt) x + dt*f(x,u);

% LQR controller
Ak = eye(n)+dt*A;
Bk = B*dt;
K = place(Ak,Bk,0.8);

% Prefilter
Kr = inv((eye(n)-Ak+Bk*K)\Bk);

eig(Ak-Bk*K)


%% Gaussian Process

% GP hyperparameters
sigmaf  = 10;        % output variance (std)
lambda  = 10;        % length scale
sigman  = sigmaw*dt; % stddev of measurement noise
maxsize = 100;       % maximum number of points in the dictionary

% create GP object
gp = GP(sigmaf, sigman, lambda, maxsize);


%% Simulate

% define input
r = @(t) 3;
r = @(t) sin(10*t);
% r = @(t) sin(10*t); %+ 2*sin(20*t) + 10*exp(-t);

nr = 1;

tf = 1;

x0 = 0;

% initialize variables to store simulation results
out.t = 0:dt:tf;
out.x = [x0 zeros(n,length(out.t)-1)];
out.xhat = [x0 zeros(n,length(out.t)-1)];
out.xnom = [x0 zeros(n,length(out.t)-1)];
out.u = zeros(m,length(out.t)-1);
out.r = zeros(nr,length(out.t)-1);


gptrue = @(x) Bd \ ( fd_true(x(1),x(2),dt,false) - fd(x(1),x(2),dt) );


for i = 1:numel(out.t)-1
    out.t(i)
    
    % read new reference
    out.r(:,i) = r(out.t(i));
    
    % calculate control input
    out.u(:,i) = Kr*out.r(:,i) - K*out.xhat(:,i);
    
    % simulate real model
    out.x(:,i+1) = fd_true(out.x(:,i),out.u(:,i),dt,true);
    
    % measure data
    out.xhat(:,i+1) = out.x(:,i+1); % perfect observer
    
    % update GP model
    out.xnom(:,i+1) = fd(out.xhat(:,i),out.u(:,i),dt);
    ddata = Bd \ (out.xhat(:,i+1) - out.xnom(:,i+1));
    gp.add( [out.xnom(:,i);out.u(:,i)], ddata );
    
    1;
    % scatter3( gp.dict.X(1,:), gp.dict.X(2,:), gp.dict.Y, 'filled' ) 
%     ddata
%     gptrue( [out.xhat(:,i),out.u(:,i)] )
%     d_true( out.xhat(:,i),out.u(:,i) )*dt
    
end


%% Evaluate results

% true GP function that is ment to be learned
% gptrue = @(x) Bd \ ( fd_true(x(1),x(2),dt,false) - fd(x(1),x(2),dt) );

% plot prediction bias and variance
gp.plot2d( gptrue, ...
           [min(out.x),max(out.x)], ...
           [min(out.u),max(out.u)] );
       
       
figure; hold on; grid on;
plot(out.t(1:end-1), out.r, 'DisplayName', 'r(t)')
plot(out.t, out.x, 'DisplayName', 'x(t)')
legend;




