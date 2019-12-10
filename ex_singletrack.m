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
tf = 50;     % simulation time


%% True Dynamics Model
%--------------------------------------------------------------------------
%   xk+1 = fd(xk,uk) + Bd * ( d(zk) + w ),    
%
%       where: zk = Bz*xk,
%              d ~ N(mean_d(zk),var_d(zk))
%              w ~ N(0,sigmaw)
%------------------------------------------------------------------

% define model (mean and variance) for true disturbance
% mu_d  = @(z) 0;
% var_d = @(z) 0;
d = @(z)deal(0,0);
sigmaw = 0;
trueModel = MotionModelGP_SingleTrack(d,sigmaw);


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
estModel = MotionModelGP_SingleTrackNominal(@d_GP.eval, sigmaw);

% create nominal dynamics model (no disturbance)
nomModel = MotionModelGP_SingleTrackNominal(@(z)deal(0,0), 0); 


%% Controller

% -------------------------------------------------------------------------
%                               (TODO)
% LQR CONTROLLER:
% -------------------------------------------------------------------------
% n = estModel.n;
% m = estModel.m;
% % % % [A,B] = estModel.linearize();
% % % % Ak = eye(n)+dt*A;
% % % % Bk = B*dt;
% % % % Ck=[0 1 0 0; 0 0 1 0; 0 0 0 1];
% % % % Q = 1e3*eye(4);
% % % % R = 1;
% % % % [~,~,K] = dare(Ak,Bk,Q,R);
% % % % % Prefilter
% % % % Kr = pinv(Ck/(eye(n)-Ak+Bk*K)*Bk);
% % % % % check eigenvalues
% % % % eig(Ak-Bk*K);



% -------------------------------------------------------------------------
% Create perception model (in this case is the saved track points)
% -------------------------------------------------------------------------
x0  = [0;0];
th0 = 0;
w = 6;
trackdata = {
     's',14;
     'c',{15,-90};
     's',5;
     'c',{4,90};
     'c',{4,-90};
     's',5;
     'c',{3.5,-180};
     'c',{3.5,180};
     'c',{3.5,-90};
     's',2;
     'c',{3.5,-120};
     's',10;
     'c',{10,120};
     's',10;
     'c',{5,90};
     's',5;
     'c',{5,150};
     's',5;
     'c',{3.2,-180};
     's',12;
     'c',{10,-150};
     's',12.3;      
     'c',{12,-90}; 
};
track = RaceTrack(trackdata, x0, th0, w);
% TEST: [Xt, Yt, PSIt, Rt] = track.getTrackInfo(1000)


trackAnim = SingleTrackAnimation(track);
trackAnim.initGraphics()



% -------------------------------------------------------------------------
% Nonlinear Model Predictive Controller
% -------------------------------------------------------------------------

% define cost function
n  = estModel.n;
m  = estModel.m;
ne = 1; % traveled track distance
N = 10; % prediction horizon
Q = diag([1000 1000 100]);
Qf= Q;
R = diag([0 0 0.1 0 -1]);
Ck = [eye(3), zeros(3,7)];

% define cost functions
fo   = @(t,mu_x,var_x,u,e,r) costFunction(mu_x, var_x, u, e, track);            % e = track distance
fend = @(t,mu_x,var_x,e,r)   costFunction(mu_x, var_x, zeros(5,1), e, track);   % end cost function

% define dynamics
f    = @(mu_x,var_x,u) estModel.xkp1(mu_x, var_x, u, dt);

% define additional constraints
h    = @(x,u,e) [];
g    = @(x,u,e) [-u(1)-deg2rad(30);  % delta >= -30deg
                  u(1)-deg2rad(30);  % delta <=  30 deg
                 -u(2)+1;            % gear  >= 1
                  u(2)-5;            % gear  <= 5
                 -u([3,4,5]);        % break, breakdist, acc >= 0
                  u([3,4,5])-1;];    % break, breakdist, acc <= 1

% Initialize NMPC object;
mpc = NMPC(f, h, g, n, m, ne, fo, fend, N, dt);
mpc.tol     = 1e-2;
mpc.maxiter = 10;

% TEST NMPC
% x0 = 10;
% t0 = 0;
% r  = @(t)2;    % desired trajectory
% u0 = mpc.optimize(x0, t0, r );



%% Simulate

% define input
% nr = size(Ck,1); % dimension of r(t)

% initial state
x0 = [10;0;0;0;0;0;0;0;0;0]; % initial value for integration

% initialize variables to store simulation results
out.t = 0:dt:tf;
out.x = [x0 zeros(n,length(out.t)-1)];
out.xhat = [x0 zeros(n,length(out.t)-1)];
out.xnom = [x0 zeros(n,length(out.t)-1)];
out.u = zeros(m,length(out.t)-1);
% out.r = zeros(nr,length(out.t)-1);


% deactivate GP evaluation in the prediction
d_GP.isActive = false;

for i = 1:numel(out.t)-1
    disp(out.t(i))
    
    % ---------------------------------------------------------------------
    % Calculate control input
    % ---------------------------------------------------------------------
    % LQR controller
    % % out.u(:,i) = Kr*out.r(:,i) - K*out.xhat(:,i);
    % NMPC controller
    [x0_opt, u_opt, e_opt] = mpc.optimize(out.xhat(:,i), out.t(i), 0);
    out.u(:,i) = u_opt(:,1);
    
    % ---------------------------------------------------------------------
    % simulate real model
    % ---------------------------------------------------------------------
    [mu_xkp1,var_xkp1] = trueModel.xkp1(out.x(:,i),zeros(trueModel.n),out.u(:,i),dt);
    out.x(:,i+1) = mvnrnd(mu_xkp1, var_xkp1, 1)';
    
    
    % ---------------------------------------------------------------------
    % measure data
    % ---------------------------------------------------------------------
    out.xhat(:,i+1) = out.x(:,i+1); % perfect observer
    
    
    % ---------------------------------------------------------------------
    % calculate nominal model
    % ---------------------------------------------------------------------
    out.xnom(:,i+1) = nomModel.xkp1(out.xhat(:,i),zeros(nomModel.n),out.u(:,i),dt);
    
    
    % ---------------------------------------------------------------------
    % add data to GP model
    % ---------------------------------------------------------------------
    if mod(i-1,1)==0
        % calculate disturbance (error between measured and nominal)
        d_est = estModel.Bd \ (out.xhat(:,i+1) - out.xnom(:,i+1));
        % select subset of coordinates that will be used in GP prediction
        zhat = estModel.Bz * out.xhat(:,i);
        % add data point to the GP dictionary
        d_GP.add(zhat,d_est);
    end
    
    if d_GP.N > 20 && out.t(i) > 3
        % d_GP.isActive = true;
    end
    
    % check if these values are the same:
    % d_est == mu_d(zhat) == ([mud,~]=trueModel.d(zhat); mud*dt)
end


%% Evaluate results
close all;



%% Help functions

function cost = costFunction(mu_x, var_x, u, dist, track)
    q_l   = 1e1; % penalization of lag error
    q_c   = 1e2; % penalization of contouring error
    q_r   = 1e4; % penalization when vehicle is outside track
    q_br  = 1e2; % penalization of breaking
    q_acc = 1e1; % reward for accelerating
    q_v   = 1e1; % reward high velocities
    q_d   = 1e5; % reward 

    % get information (x,y,track radius and track orientation) of the point 
    % in the track that corresponds to a traveled distance of 'dist' meters.
    track_dist = mu_x(11);
    [pos_c, psi_c, R_c] = track.getTrackInfo(track_dist);
    
    % ---------------------------------------------------------------------
    % cost of contour and lag error
    % ---------------------------------------------------------------------
    % rotation to a frame with x-axis tangencial to the track (T frame)
    A_TI = [cos(psi_c) -sin(psi_c);    
            sin(psi_c)  cos(psi_c)];
    % error in the inertial coordinates
    I_error = pos_c - mu_x(1:2);       
    % error in the T frame [lag_error; contouring_error]
    T_error = A_TI * I_error;          
    
    cost_contour = q_c*T_error(2)^2;
    cost_lag     = q_l*T_error(1)^2;
    
    % ---------------------------------------------------------------------
    % cost for being outside track
    % ---------------------------------------------------------------------
    % is the vehicle outside the track?
    isOusideTrack = (norm(I_error) - R_c) > 0;
    cost_outside = isOusideTrack*q_r*(norm(I_error)-R_c)^2;
    
    % ---------------------------------------------------------------------
    % reward high velocities
    % ---------------------------------------------------------------------
    cost_vel = - q_v*mu_x(3)^2;
    
    % ---------------------------------------------------------------------
    % cost for high inputs - penalize high inputs
    % ---------------------------------------------------------------------
    brake_force = u(3);
    acc_pedal   = u(5);
    cost_inputs = q_br*brake_force^2 - q_acc*acc_pedal^2;
    
    % ---------------------------------------------------------------------
    % reward track velocities
    % ---------------------------------------------------------------------
    track_vel = u(6);
    cost_dist = -q_d*track_vel;
    
    % ---------------------------------------------------------------------
    % Calculate final cost
    % ---------------------------------------------------------------------
    cost = cost_contour + cost_lag + cost_outside + cost_inputs + cost_vel + cost_dist;
end


