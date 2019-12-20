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

dt = 0.1;   % simulation timestep size
tf = 50;     % simulation time


%% True Dynamics Model
%--------------------------------------------------------------------------
%   xk+1 = fd(xk,uk) + Bd * ( d(zk) + w ),    
%
%       where: zk = Bz*xk,
%              d ~ N(mean_d(zk),var_d(zk))
%              w ~ N(0,var_w)
%------------------------------------------------------------------

% define model (mean and variance) for true disturbance
% mu_d  = @(z) 0;
% var_d = @(z) 0;
d = @(z)deal(0,0);
var_w = (1/3)^2;    %diag([(1/3)^2 (1/3)^2 (deg2rad(1)/3)^2]);
trueModel = MotionModelGP_SingleTrackNominal(d,var_w);


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
estModel = MotionModelGP_SingleTrackNominal(@d_GP.eval, var_w);

% create nominal dynamics model (no disturbance)
nomModel = MotionModelGP_SingleTrackNominal(@(z)deal(0,0), 0*var_w); 


%% Controller

% -------------------------------------------------------------------------
%                               (TODO)
% LQR CONTROLLER:
% -------------------------------------------------------------------------
% % % % [A,B] = estModel.linearize();
% % % % Ak = eye(estModel.n)+dt*A;
% % % % Bk = B*dt;
% % % % Ck=[0 1 0 0; 0 0 1 0; 0 0 0 1];
% % % % Q = 1e3*eye(estModel.n);
% % % % R = 1;
% % % % [~,~,K] = dare(Ak,Bk,Q,R);
% % % % % Prefilter
% % % % Kr = pinv(Ck/(eye(estModel.n)-Ak+Bk*K)*Bk);
% % % % % check eigenvalues
% % % % eig(Ak-Bk*K);



% -------------------------------------------------------------------------
% Create perception model (in this case is the saved track points)
% -------------------------------------------------------------------------
[trackdata, x0, th0, w] = RaceTrack.loadTrack_01();
track = RaceTrack(trackdata, x0, th0, w);
% TEST: [Xt, Yt, PSIt, Rt] = track.getTrackInfo(1000)
%       trackAnim = SingleTrackAnimation(track,mpc.N);
%       trackAnim.initGraphics()

% -------------------------------------------------------------------------
% Nonlinear Model Predictive Controller
% -------------------------------------------------------------------------

% define cost function
n  = estModel.n;
m  = estModel.m;
ne = 0;

N = 10; % prediction horizon

lookahead = dt*N

% define cost functions
fo   = @(t,mu_x,var_x,u,e,r) costFunction(mu_x, var_x, u, track);            % e = track distance
fend = @(t,mu_x,var_x,e,r)   10 * costFunction(mu_x, var_x, zeros(m,1), track);   % end cost function

% define dynamics
f  = @(mu_x,var_x,u) estModel.xkp1(mu_x, var_x, u, dt);
% define additional constraints
h  = @(x,u,e) [];
g  = @(x,u,e) [];
u_lb = [-deg2rad(30);  % delta >= -10deg
         -1;           % wheel torque gain >= -1
         5];           % track velocity >= 0
u_ub = [deg2rad(30);   % delta <=  10 deg
        1;             % wheel torque gain <= 1
        30];           % track velocity <= 1

% Initialize NMPC object;
mpc = NMPC(f, h, g, u_lb, u_ub, n, m, ne, fo, fend, N, dt);
mpc.tol     = 1e-2;
mpc.maxiter = 30;

% TEST NMPC
% x0 = 10;
% t0 = 0;
% r  = @(t)2;    % desired trajectory
% u0 = mpc.optimize(x0, t0, r );



%% Prepare simulation

% ---------------------------------------------------------------------
% Prepare simulation (initialize vectors, initial conditions and setup
% animation
% ---------------------------------------------------------------------

% define variable sizes
true_n = trueModel.n;
true_m = trueModel.m;
est_n = estModel.n;
est_m = estModel.m;

% initial state
x0 = [10;0;0; 10;0;0; 0;];   % true initial state
x0(end) = track.getTrackDistance(x0(1:2)); % get initial track traveled distance

% change initial guess for mpc solver. Set initial track velocity as
% initial vehicle velocity (this improves convergence speed a lot)
mpc.uguess(end,:) = x0(4)*2;

% define simulation time
out.t = 0:dt:tf;            % time vector
kmax = length(out.t)-1;     % steps to simulate

% initialize variables to store simulation results
out.x              = [x0 NaN(true_n,kmax)];
out.xhat           = [x0 NaN(est_n, kmax)];
out.xnom           = [x0 NaN(est_n, kmax)];
out.u              =     NaN(est_m, kmax);
out.x_ref          = NaN(2,     mpc.N+1, kmax);
out.mu_x_pred_opt  = NaN(mpc.n, mpc.N+1, kmax);
out.var_x_pred_opt = NaN(mpc.n, mpc.n, mpc.N+1, kmax);
out.u_pred_opt     = NaN(mpc.m, mpc.N,   kmax);


% start animation
trackAnim = SingleTrackAnimation(track,out.mu_x_pred_opt,out.u_pred_opt,out.x_ref);
trackAnim.initTrackAnimation();
trackAnim.initScope();

% deactivate GP evaluation in the prediction
d_GP.isActive = false;

%% Start simulation
% ---------------------------------------------------------------------
% Start simulation
% ---------------------------------------------------------------------
ki = 1;
% mpc.uguess = out.u_pred_opt(:,:,ki);

for k = ki:kmax
    disp(out.t(k))
    
    % ---------------------------------------------------------------------
    % LQR controller
    % ---------------------------------------------------------------------
    % % out.u(:,i) = Kr*out.r(:,i) - K*out.xhat(:,i);
    
    % ---------------------------------------------------------------------
    % NPMC controller
    % ---------------------------------------------------------------------
    % calculate optimal input
    [x0_opt, u_opt, e_opt] = mpc.optimize(out.xhat(:,k), out.t(k), 0);
    out.u(:,k) = u_opt(:,1);
    sprintf('\nSteering angle: %d\nTorque gain: %.1f\nTrack vel: %.1f\n',rad2deg(out.u(1,k)),out.u(2,k),out.u(3,k))

    % ---------------------------------------------------------------------
    % store data and plot
    % ---------------------------------------------------------------------
    % get optimal state predictions from optimal input and current state
    out.u_pred_opt(:,:,k) = u_opt;
    [out.mu_x_pred_opt(:,:,k),out.var_x_pred_opt(:,:,:,k)] = mpc.predictStateSequence(out.xhat(:,k), zeros(estModel.n), u_opt);
    % get target track distances from predictions (last state)
    out.x_ref(:,:,k) = track.getTrackInfo(out.mu_x_pred_opt(end,:,k));
    
    % update track animation
    trackAnim.mu_x_pred_opt = out.mu_x_pred_opt;
    trackAnim.u_pred_opt = out.u_pred_opt;
    trackAnim.x_ref      = out.x_ref;
    trackAnim.updateTrackAnimation(k);
    trackAnim.updateScope(k);
    
    % ---------------------------------------------------------------------
    % simulate real model
    % ---------------------------------------------------------------------
    [mu_xkp1,var_xkp1] = trueModel.xkp1(out.x(:,k),zeros(trueModel.n),out.u(:,k),dt);
    out.x(:,k+1) = mvnrnd(mu_xkp1, var_xkp1, 1)';
    
    
    % ---------------------------------------------------------------------
    % measure data
    % ---------------------------------------------------------------------
    out.xhat(:,k+1) = out.x(:,k+1); % perfect observer
    % get traveled distance, given vehicle coordinates (this is the 11th
    % state of the nominal model)
    out.xhat(end,k+1) = track.getTrackDistance(out.xhat([1,2],k+1));
    
    
    % ---------------------------------------------------------------------
    % Safety
    % ---------------------------------------------------------------------
    V_vx = out.xhat(4,k+1);
    V_vy = out.xhat(5,k+1);
    beta = atan2(V_vy,V_vx);
    if V_vx < 0
        error('Vehicle is driving backwards... aborting');
    end
    if abs(rad2deg(beta)) > 80
        error('Vehicle has a huge sideslip angle... aborting')
    end    
    
    % ---------------------------------------------------------------------
    % calculate nominal model
    % ---------------------------------------------------------------------
    out.xnom(:,k+1) = nomModel.xkp1(out.xhat(:,k),zeros(nomModel.n),out.u(:,k),dt);
    
    
    % ---------------------------------------------------------------------
    % add data to GP model
    % ---------------------------------------------------------------------
    % if mod(k-1,1)==0
    %     % calculate disturbance (error between measured and nominal)
    %     d_est = estModel.Bd \ (out.xhat(:,k+1) - out.xnom(:,k+1));
    %     % select subset of coordinates that will be used in GP prediction
    %     zhat = estModel.Bz * out.xhat(:,k);
    %     % add data point to the GP dictionary
    %     d_GP.add(zhat,d_est);
    %     d_GP.updateModel();
    % end
    % 
    % if d_GP.N > 50 && out.t(k) > 3
    %     % d_GP.isActive = true;
    % end
    
    % check if these values are the same:
    % d_est == mu_d(zhat) == ([mud,~]=trueModel.d(zhat); mud*dt)
   
end


%% Show animation
close all;

% start animation
trackAnim = SingleTrackAnimation(track,out.mu_x_pred_opt,out.u_pred_opt,out.x_ref);
trackAnim.initTrackAnimation();
% trackAnim.initScope();
for k = 1:kmax
    trackAnim.updateTrackAnimation(k);
    % trackAnim.updateScope(k);
    pause(0.1);
end

%% Record video

trackAnim.recordvideo(fullfile('simresults','trackAnimVideo'),'Motion JPEG AVI');


%% Help functions

function cost = costFunction(mu_x, var_x, u, track)

    % Track oriented penalization
    q_l   = 1e2; % penalization of lag error
    q_c   = 1e2; % penalization of contouring error
    q_o   = 1e1; % penalization for orientation error
    q_d   = 1e0; % reward high track centerline velocites
    q_r   = 1e3; % penalization when vehicle is outside track
    
    % state and input penalization
    q_st  = 1*1e1; % penalization of steering
    q_br  = 0*1e0; % penalization of breaking
    q_acc = 0*1e0; % reward for acceleration
    q_v   = 1*1e0; % reward high absolute velocities

    % label inputs and outputs
    I_x        = mu_x(1);  % x position in global coordinates
    I_y        = mu_x(2);  % y position in global coordinates
    psi        = mu_x(3);  % yaw
    V_vx       = mu_x(4);  % x velocity in vehicle coordinates
    V_vy       = mu_x(5);  % x velocity in vehicle coordinates
    track_dist = mu_x(7);  % track velocity
    delta      = u(1);     % steering angle rad2deg(delta)
    T          = u(2);     % torque gain (1=max.acc, -1=max.braking)
    track_vel  = u(3);     % track velocity
    

    % ---------------------------------------------------------------------
    % cost of contour, lag and orientation error
    % ---------------------------------------------------------------------

    % get lag, contour, offroad and orientation error of the vehicle w.r.t.
    % a point in the trajectory that is 'track_dist' far away from the 
    % origin along the track centerline (traveled distance)
    [lag_error, countour_error, offroad_error, orientation_error] = ...
        track.getVehicleDeviation([I_x;I_y], psi, track_dist);
    
    cost_contour     = q_c * countour_error^2;
    cost_lag         = q_l * lag_error^2;
    cost_orientation = q_o * orientation_error^2;
    
    % ---------------------------------------------------------------------
    % cost for being outside track
    % ---------------------------------------------------------------------
    % is the vehicle outside the track?
    cost_outside = q_r * offroad_error^2;
    
    % ---------------------------------------------------------------------
    % reward high velocities only if inside track
    % ---------------------------------------------------------------------
    cost_vel = -q_v * norm([V_vx; V_vy]);
    
    % ---------------------------------------------------------------------
    % reward high track velocities
    % ---------------------------------------------------------------------
    cost_dist = -q_d * track_vel;
    
    % ---------------------------------------------------------------------
    % reward acceleration and penalize braking and steering
    % ---------------------------------------------------------------------
    cost_inputs = - (T>0)*q_acc*T^2 + (T<0)*q_br*(T)^2 + q_st*(delta)^2 ;
    
    % ---------------------------------------------------------------------
    % Calculate final cost
    % ---------------------------------------------------------------------
    cost = cost_contour + cost_lag + cost_orientation + cost_dist + cost_outside + cost_inputs + cost_vel;
    
% fprintf('Contribution to cost:\n')
% fprintf('   cost_contour:%.1f\n',cost_contour/cost*100);
% fprintf('   cost_lag:%.1f\n',cost_lag/cost*100);
% fprintf('   cost_orientation:%.1f\n',cost_orientation/cost*100);
% fprintf('   cost_dist:%.1f\n',cost_dist/cost*100);
% fprintf('   cost_outside:%.1f\n',cost_outside/cost*100);
% fprintf('   cost_inputs:%.1f\n',cost_inputs/cost*100);
% fprintf('   cost_vel:%.1f\n',cost_vel/cost*100);
end
