%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -

%   Control of a Race Car in a Race Track using Gaussian Process Model Predictive Control:
%------------------------------------------------------------------

clear all; close all; clc;


%--------------------------------------------------------------------------
% Quick Access Simulation and controller parameters
%------------------------------------------------------------------
dt = 0.15;       % simulation timestep size
tf = 15*12;       % simulation time
maxiter = 30;   % max NMPC iterations per time step
N = 10;         % NMPC prediction horizon

loadPreTrainedGP = false;
GPfile = fullfile(pwd,'/simresults/20-01-15-out-GP-with-GP-optimized.mat');
useGP = false;
trainGPonline = true;
useParallel = true;


% display info
lookahead = dt*N;
fprintf('\nPrediction lookahead: %.1f [s]\n',lookahead);



%% Create True Dynamics Simulation Model
%--------------------------------------------------------------------------
%   xk+1 = fd_true(xk,uk) + Bd * ( w ),    
%
%       where: w ~ N(0,var_w)
%------------------------------------------------------------------

% define noise for true disturbance
var_w = diag([(1/3)^2 (1/3)^2 (deg2rad(3)/3)^2]);
% var_w = zeros(3);

% create true dynamics model
trueModel = MotionModelGP_SingleTrack_true( [], var_w);
% trueModel = MotionModelGP_SingleTrack_nominal(d,var_w);


%% Create Estimation Model and Nominal Model

% -------------------------------------------------------------------------
%  Create nominal model (no disturbance):  
%       xk+1 = fd_nom(xk,uk)
% -------------------------------------------------------------------------

nomModel = MotionModelGP_SingleTrack_nominal( [], [] ); 
% nomModel = MotionModelGP_SingleTrack_true( [], [] );

nomModel.analyseSingleTrack();


% -------------------------------------------------------------------------
%  Create adaptive dynamics model 
%  (unmodeled dynamics will be estimated by Gaussian Process GP)
%       xk+1 = fd_nom(xk,uk) + Bd * ( d_GP(zk) + w )
% -------------------------------------------------------------------------

if ~loadPreTrainedGP
    % GP input dimension
    gp_n = MotionModelGP_SingleTrack_nominal.nz;
    % GP output dimension
    gp_p = MotionModelGP_SingleTrack_nominal.nd;

    % GP hyperparameters
    var_f   = repmat(0.01,[gp_p,1]);    % output variance
    var_n   = diag(var_w/3);              % measurement noise variance
    M       = repmat(diag([1e0,1e0,1e0,1e0,1e0].^2),[1,1,gp_p]);     % length scale
    maxsize = 300; % maximum number of points in the dictionary

    % create GP object
    d_GP = GP(gp_n, gp_p, var_f, var_n, M, maxsize);
else
    load(fullfile(pwd,'/simresults/20-01-15-out-GP-WORKED-optimized.mat')); %,'d_GP'
    fprintf('\nGP model loaded succesfuly\n\n')
end

% create nominal model with GP model as d(zk)
estModel = MotionModelGP_SingleTrack_nominal(@d_GP.eval, var_w);
% estModel = MotionModelGP_SingleTrack_true(@d_GP.eval, var_w);


%% Initialize Controller

% -------------------------------------------------------------------------
% Create perception model (in this case is the saved track points)
% this is needed to for the MPC cost function
% -------------------------------------------------------------------------
[trackdata, x0, th0, w] = RaceTrack.loadTrack_02();
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

% define cost functions
fo   = @(t,mu_x,var_x,u,e,r) costFunction(mu_x, var_x, u, track);            % e = track distance
fend = @(t,mu_x,var_x,e,r)   2 * costFunction(mu_x, var_x, zeros(m,1), track);   % end cost function

% define dynamics
 f  = @(mu_x,var_x,u) estModel.xkp1(mu_x, var_x, u, dt);
%f  = @(mu_x,var_x,u) trueModel.xkp1(mu_x, var_x, u, dt);
% define additional constraints
h  = @(x,u,e) [];
g  = @(x,u,e) [];
u_lb = [-deg2rad(20);  % >= steering angle
         -1;           % >= wheel torque gain
         5];           % >= centerline track velocity
u_ub = [deg2rad(20);   % <= steering angle
        1;             % <= wheel torque gain
        30];           % <= centerline track velocity 

% Initialize NMPC object;
mpc = NMPC(f, h, g, u_lb, u_ub, n, m, ne, fo, fend, N, dt);
mpc.tol     = 1e-2;
mpc.maxiter = maxiter;



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
x0 = [5;0;0; 10;0;0; 0];   % true initial state
x0(end) = track.getTrackDistance(x0(1:2)); % get initial track traveled distance

% change initial guess for mpc solver. Set initial track velocity as
% initial vehicle velocity (this improves convergence speed a lot)
mpc.uguess(end,:) = x0(4)*2;

% define simulation time
out.t = 0:dt:tf;            % time vector
kmax = length(out.t)-1;     % steps to simulate

% initialize variables to store simulation results
out.x              = [x0 NaN(true_n,kmax)];             % true states
out.xhat           = [x0 NaN(est_n, kmax)];             % state estimation
out.xnom           = [x0 NaN(est_n, kmax)];             % predicted nominal state
out.u              =     NaN(est_m, kmax);              % applied input
out.x_ref          = NaN(2,     mpc.N+1, kmax);         % optimized reference trajectory
out.mu_x_pred_opt  = NaN(mpc.n, mpc.N+1, kmax);         % mean of optimal state prediction sequence
out.var_x_pred_opt = NaN(mpc.n, mpc.n, mpc.N+1, kmax);  % variance of optimal state prediction sequence
out.u_pred_opt     = NaN(mpc.m, mpc.N,   kmax);         % open-loop optimal input prediction


% start animation
trackAnim = SingleTrackAnimation(track, out.mu_x_pred_opt, out.var_x_pred_opt, out.u_pred_opt, out.x_ref);
trackAnim.initTrackAnimation();
trackAnim.initScope();
drawnow;

% deactivate GP evaluation in the prediction
d_GP.isActive = useGP;
fprintf('\nGP active? %s\n\n',string(useGP))




%% Start simulation

ki = 1;
% ki = 310;
% mpc.uguess = out.u_pred_opt(:,:,ki);


for k = ki:kmax
    disp('------------------------------------------------------')
    fprintf('time: %.3f [s]\n',out.t(k))
    
    % ---------------------------------------------------------------------
    % LQR controller
    % ---------------------------------------------------------------------
    % % out.u(:,i) = Kr*out.r(:,i) - K*out.xhat(:,i);
    
    % ---------------------------------------------------------------------
    % NPMC controller
    % ---------------------------------------------------------------------
    % calculate optimal input
    [u_opt, e_opt] = mpc.optimize(out.xhat(:,k), out.t(k), 0, useParallel);
    out.u(:,k) = u_opt(:,1);
    sprintf('\nSteering angle: %d\nTorque gain: %.1f\nTrack vel: %.1f\n',rad2deg(out.u(1,k)),out.u(2,k),out.u(3,k))

    % ---------------------------------------------------------------------
    % Calculate predicted trajectory from optimal open-loop input sequence 
    % and calculate optimized reference trajectory for each prediction
    % ---------------------------------------------------------------------
    % get optimal state predictions from optimal input and current state
    out.u_pred_opt(:,:,k) = u_opt;
    [out.mu_x_pred_opt(:,:,k),out.var_x_pred_opt(:,:,:,k)] = mpc.predictStateSequence(out.xhat(:,k), zeros(estModel.n), u_opt);
    % get target track distances from predictions (last state)
    out.x_ref(:,:,k) = track.getTrackInfo(out.mu_x_pred_opt(end,:,k));
    
    % ---------------------------------------------------------------------
    % update race animation and scopes
    % ---------------------------------------------------------------------
    trackAnim.mu_x_pred_opt  = out.mu_x_pred_opt;
    trackAnim.var_x_pred_opt = out.var_x_pred_opt;
    trackAnim.u_pred_opt     = out.u_pred_opt;
    trackAnim.x_ref          = out.x_ref;
    trackAnim.updateTrackAnimation(k);
    trackAnim.updateScope(k);
    drawnow;
    
    % ---------------------------------------------------------------------
    % Simulate real model
    % ---------------------------------------------------------------------
    [mu_xkp1,var_xkp1] = trueModel.xkp1(out.x(:,k),zeros(trueModel.n),out.u(:,k),dt);
    % out.x(:,k+1) = mvnrnd(mu_xkp1, var_xkp1, 1)';
    out.x(:,k+1) = mu_xkp1;
    
    % ---------------------------------------------------------------------
    % Measure data
    % ---------------------------------------------------------------------
    out.xhat(:,k+1) = out.x(:,k+1); % perfect observer
    % get traveled distance, given vehicle coordinates
    out.xhat(end,k+1) = track.getTrackDistance( out.xhat([1,2],k+1) , out.xhat(end,k) );
    
    
    % ---------------------------------------------------------------------
    % Lap timer
    % ---------------------------------------------------------------------
    [laptimes, idxnewlaps] = RaceTrack.getLapTimes(out.xhat(end,:),dt);
    if any(k==idxnewlaps)
        RaceTrack.dispLapTimes(laptimes);
    end
    
    
    % ---------------------------------------------------------------------
    % Safety - Stop simulation in case vehicle is completely unstable
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
    % Add data to GP model
    % ---------------------------------------------------------------------
    if mod(k-1,1)==0
        % calculate disturbance (error between measured and nominal)
        d_est = estModel.Bd \ (out.xhat(:,k+1) - out.xnom(:,k+1));
        % d_est = estModel.Bd \ (mu_xkp1 - out.xnom(:,k+1));
        % select subset of coordinates that will be used in GP prediction
        zhat = [ estModel.Bz_x * out.xhat(:,k); estModel.Bz_u * out.u(:,k) ];
        % add data point to the GP dictionary
        if trainGPonline
            d_GP.add(zhat,d_est');
            d_GP.updateModel();
        end
        
        fprintf('Prediction Error norm WITHOUT GP: %f\n',norm(d_est));
        disp(d_est)
        fprintf('Prediction Error norm WITH    GP: %f\n',norm(d_est-d_GP.eval(zhat,true)));
        % fprintf('Prediction Error norm WITH    GP: %f\n',norm(d_est-estModel.d(zhat,true)));
        disp(d_est-d_GP.eval(zhat,true))
    end
    
    % if length(laptimes) >= 6
        % d_GP.isActive = true;
        % mpc.maxiter = 30;
    % end

end


% Display Lap times

[laptimes, idxnewlaps] = RaceTrack.getLapTimes(out.xhat(end,:),dt);
RaceTrack.dispLapTimes(laptimes)




return
% STOP here. Next sections are intended to be executed separately



%% Readd simulation data to GP, uddate model and optimize parameters

k = find(~isnan(out.xhat(1,:)), 1, 'last' ) - 20;

% create new instance of GP class
d_GP = GP(gp_n, gp_p, var_f, var_n, M, maxsize);

% readd points
d_est = estModel.Bd \ (out.xhat(:,2:k) - out.xnom(:,2:k));
zhat  = estModel.z( out.xhat(:,1:k-1), out.u(:,1:k-1) );
d_GP.add(zhat,d_est');

% update and optimize model
d_GP.updateModel();
d_GP.optimizeHyperParams('fmincon')

d_GP.M
d_GP.var_f
d_GP.var_n


%% Analyse learning
% ---------------------------------------------------------------------
% Check how the GP reduces the prediction error
% ---------------------------------------------------------------------

% d_GP.optimizeHyperParams('fmincon')
% d_GP.optimizeHyperParams('ga')


k = find(~isnan(out.xhat(1,:)), 1, 'last' ) - 20;

% prediction error without GP
% predErrorNOgp = estModel.Bd\(out.xhat - out.xnom);
predErrorNOgp = estModel.Bd\(out.xhat(:,1:k-1) - out.xnom(:,1:k-1));


% prediction error with trained GP
zhat  = estModel.z( out.xhat(:,1:k-1), out.u(:,1:k-1) );
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


% ---------------------------------------------------------------------
% Check in which region of the tyre dynamics we are working
% ---------------------------------------------------------------------

% % % % simulation
% % % 
% trueModel.testTyres
% 
% l_f  = 0.9;
% l_r  = 1.5;
% V_vx = out.xhat(4,:);
% V_vy = out.xhat(5,:);
% psi_dot = out.xhat(6,:);
% delta = out.u(1,:);
% a_r = atan2(V_vy-l_r.*psi_dot,V_vx);
% a_f = atan2(V_vy+l_f.*psi_dot,V_vx) - [delta 0];
% 
% figure('Color','w'); hold on; grid on;
% plot(rad2deg(a_r))
% plot(rad2deg(a_f))
% ylabel('slip angle')
% xlabel('time step')



%% Show animation
close all;

% start animation
trackAnim = SingleTrackAnimation(track, out.mu_x_pred_opt, out.var_x_pred_opt, out.u_pred_opt, out.x_ref);
trackAnim.initTrackAnimation();
% trackAnim.initScope();
for k=1:kmax
    if ~ trackAnim.updateTrackAnimation(k)
        break;
    end
    % trackAnim.updateScope(k);
%     pause(0.15);
    drawnow;
end


%% Record video

FrameRate = 7;
videoName = fullfile('simresults',sprintf('trackAnimVideo-%s',date));
videoFormat = 'Motion JPEG AVI';
trackAnim.recordvideo(videoName, videoFormat, FrameRate);


%% Cost function for the MPC controller

function cost = costFunction(mu_x, var_x, u, track)

    % Track oriented penalization
    q_l   = 50;     % penalization of lag error
    q_c   = 20;     % penalization of contouring error
    q_o   = 5;      % penalization for orientation error
    q_d   = -3;     % reward high track centerline velocites
    q_r   = 100;    % penalization when vehicle is outside track
    
    % state and input penalization
    q_v      = -0;  % reward high absolute velocities
    q_st     =  0;  % penalization of steering
    q_br     =  0;  % penalization of breaking
    q_psidot =  8;  % penalize high yaw rates
    q_acc    = -0;  % reward for accelerating

    % label inputs and outputs
    I_x        = mu_x(1);  % x position in global coordinates
    I_y        = mu_x(2);  % y position in global coordinates
    psi        = mu_x(3);  % yaw
    V_vx       = mu_x(4);  % x velocity in vehicle coordinates
    V_vy       = mu_x(5);  % x velocity in vehicle coordinates
    psidot     = mu_x(6);
    track_dist = mu_x(7);  % track centerline distance
    delta      = u(1);     % steering angle rad2deg(delta)
    T          = u(2);     % torque gain (1=max.acc, -1=max.braking)
    track_vel  = u(3);     % track centerline velocity
    

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
    % % apply smooth barrier function (we want: offroad_error < 0). 
    % alpha = 40; % smoothing factor... the smaller the smoother
    % offroad_error = (1+exp(-alpha*(offroad_error+0.05))).^-1;
    gamma = 1000;
    lambda = -0.1;
    offroad_error = 5*(sqrt((4+gamma*(lambda-offroad_error).^2)/gamma) - (lambda-offroad_error));

    % % CHECK SMOOTH TRANSITION
    % x = -0.5:0.01:0.5
    % % Smooth >=0 boolean function
    % alpha = 40; % the larger the sharper the clip function
    % y = (1+exp(-alpha*(x+0.05))).^-1 + exp(x);
    % gamma = 10000;
    % lambda = -0.2;
    % y = 0.5*(sqrt((4+gamma*(lambda-x).^2)/gamma) - (lambda-x));
    % figure; hold on; grid on;
    % plot(x,y)
    cost_outside = q_r * offroad_error^2;
    
    % ---------------------------------------------------------------------
    % reward high velocities
    % ---------------------------------------------------------------------
    cost_vel = q_v * norm([V_vx; V_vy]);
    
    % ---------------------------------------------------------------------
    % penalize high yaw rates
    % ---------------------------------------------------------------------
    cost_psidot = q_psidot * psidot^2;
    
    % ---------------------------------------------------------------------
    % reward high track velocities
    % ---------------------------------------------------------------------
    cost_dist = q_d * track_vel;
    
    % ---------------------------------------------------------------------
    % penalize acceleration, braking and steering
    % ---------------------------------------------------------------------
    cost_inputs = (T>0)*q_acc*T^2 + (T<0)*q_br*T^2 + q_st*(delta)^2 ;
    
    % ---------------------------------------------------------------------
    % Calculate final cost
    % ---------------------------------------------------------------------
    cost = cost_contour + ...
           cost_lag + ...
           cost_orientation + ...
           cost_dist + ...
           cost_outside + ...
           cost_inputs + ...
           cost_vel + ...
           cost_psidot;
end