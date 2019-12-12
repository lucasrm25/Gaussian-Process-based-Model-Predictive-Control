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
trueModel = MotionModelGP_SingleTrackNominal(d,sigmaw);


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
[trackdata, x0, th0, w] = RaceTrack.loadTrack_01();
track = RaceTrack(trackdata, x0, th0, w);
% TEST: [Xt, Yt, PSIt, Rt] = track.getTrackInfo(1000)


% -------------------------------------------------------------------------
% Nonlinear Model Predictive Controller
% -------------------------------------------------------------------------

% define cost function
n  = estModel.n;
m  = estModel.m;
ne = 0;

N = 8; % prediction horizon

% define cost functions
fo   = @(t,mu_x,var_x,u,e,r) costFunction(mu_x, var_x, u, track);            % e = track distance
fend = @(t,mu_x,var_x,e,r)   costFunction(mu_x, var_x, zeros(m,1), track);   % end cost function

% define dynamics
f  = @(mu_x,var_x,u) estModel.xkp1(mu_x, var_x, u, dt);
% define additional constraints
h  = @(x,u,e) [];
g  = @(x,u,e) [];
u_lb = [-deg2rad(20);  % delta_dot >= -10deg/s
         -1;           % wheel torque gain >= -1
         0;            % wheel torque distribution >= 0
         3];           % track velocity >= 0
u_ub = [deg2rad(20);   % delta_dot <=  10 deg/s
        1;             % wheel torque gain <= 1
        1;             % wheel torque distribution <= 1
        10];           % track velocity <= 1

% Initialize NMPC object;
mpc = NMPC(f, h, g, u_lb, u_ub, n, m, ne, fo, fend, N, dt);
mpc.tol     = 1e-3;
mpc.maxiter = 10;

% TEST NMPC
% x0 = 10;
% t0 = 0;
% r  = @(t)2;    % desired trajectory
% u0 = mpc.optimize(x0, t0, r );



%% Simulate

% define variable sizes
true_n = trueModel.n;
true_m = trueModel.m;
est_n = estModel.n;
est_m = estModel.m;

% initial state
true_x0 = [10;0;0; 5;0;0; 0; 0];   % true initial state
true_x0(end) = track.getTrackDistance(true_x0(1:2)); % get initial track traveled distance
est_x0  = [10;0;0; 5;0;0; 0; 0];   % initial state prior
est_x0(end) = track.getTrackDistance(est_x0(1:2)); % get initial track traveled distance

% define simulation time
out.t = 0:dt:tf;            % time vector
kmax = length(out.t)-1;     % steps to simulate

% initialize variables to store simulation results
out.x    = [true_x0 NaN(true_n,kmax)];
out.xhat = [est_x0  NaN(est_n,kmax)];
out.xnom = [est_x0  NaN(est_n,kmax)];
out.u    = NaN(est_m,kmax);
out.ref     = NaN(2,mpc.N+1,kmax);
out.estPred = NaN(2,mpc.N+1,kmax);

% deactivate GP evaluation in the prediction
d_GP.isActive = false;

% start animation
trackAnim = SingleTrackAnimation(track,mpc.N);
trackAnim.initGraphics()

plotscope = true;
if plotscope
    scopex = figure('Position',[-1006 86 957 808]);
    scopeu = figure('Position',[-1879 93 867 795]);
    plotScope(scopex,scopeu,out.xhat,out.u);
end
    
% ---------------------------------------------------------------------
% Start simulation
% ---------------------------------------------------------------------
for k = 1:kmax
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
    out.u(:,k) = u_opt(:,1); % get first input and discard last virtual input (track velocity)
    sprintf('\nSteering angle velocity: %d\nTorque gain: %.1f\nTorque dist: %.1f\nTrack vel: %.1f\n',rad2deg(out.u(1,k)),out.u(2,k),out.u(3,k),out.u(4,k))

    % ---------------------------------------------------------------------
    % plot
    % ---------------------------------------------------------------------
    estPred = mpc.predictStateSequence(out.xhat(:,k), zeros(estModel.n), u_opt);
    out.ref(:,:,k)     = track.getTrackInfo(estPred(end,:));
    out.estPred(:,:,k) = estPred(1:2,:);     % we only need X,Y
    trackAnim.estPred  = out.estPred(:,:,k);
    trackAnim.ref      = out.ref(:,:,k);
    trackAnim.updateGraphics();
    if plotscope
        refreshdata(scopex);
        refreshdata(scopeu);
        drawnow;
    end
    
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
trackAnim = SingleTrackAnimation(track,mpc.N);
trackAnim.initGraphics()
% start scope
scopex = figure('Position',[-1006 86 957 808]);
scopeu = figure('Position',[-1879 93 867 795]);
plotScope(scopex,scopeu,out.xhat,out.u);
for k = 1:kmax
    trackAnim.estPred  = out.estPred(:,:,k);
    trackAnim.ref      = out.ref(:,:,k);
    trackAnim.updateGraphics();
    drawnow
    pause(0.2);
end


%% Help functions

function plotScope(figx, figu, x, u)
%   x = [I_x              (x position in global coordinates), 
%        I_y              (y position in global coordinates), 
%        psi              (yaw angle),
%        V_vx             (longitudinal velocity in vehicle coordinates)             
%        V_vy             (lateral velocity in vehicle coordinates)
%        psi_dot          (yaw rate),
%        delta            (steering angle)
%        track_dist       (distance traveled in the track centerline)
%        ]
%
%   u = [delta_dot      (steering angle velocity), 
%        G              (gear), 
%        F_b            (brake force),
%        zeta           (brake force distribution), 
%        phi            (acc pedal position),
%        track_vel      (velocity in the track centerline)
%       ] 
    figure(figx);
    names = {'I-x','I-y','psi','V-vx','V-vy','psidot','delta','track_dist'};
    angles = [0 0 1 0 0 1 1 0];
    for i=1:numel(names)
        subplot(4,2,i);
        if angles(i)
            p = plot(rad2deg(x(i,:)),'DisplayName',names{i});
            p.XDataSource = sprintf('out.t');
            p.YDataSource = sprintf('rad2deg(out.x(%d,:))',i);
        else
            p = plot(x(i,:),'DisplayName',names{i});
            p.XDataSource = sprintf('out.t');
            p.YDataSource = sprintf('out.x(%d,:)',i);
        end
        legend('Location', 'southwest');
    end
    figure(figu);
    names = {'deltadot','T','T_dist','Track vel'};
    angles = [1 0 0 0];
    for i=1:numel(names)
        subplot(2,2,i);
        if angles(i)
            p = plot(rad2deg(u(i,:)),'DisplayName',names{i});
            p.XDataSource = sprintf('out.t(1:end-1)');
            p.YDataSource = sprintf('rad2deg(out.u(%d,:))',i);
        else
            p = plot(u(i,:),'DisplayName',names{i});
            p.XDataSource = sprintf('out.t(1:end-1)');
            p.YDataSource = sprintf('out.u(%d,:)',i);
        end
        legend('Location', 'southwest');
    end
end

function cost = costFunction(mu_x, var_x, u, track)

    % Track oriented penalization
    q_l   = 1e5; % penalization of lag error
    q_c   = 1e5; % penalization of contouring error
    q_o   = 1e2; % penalization for orientation error
    q_d   = 1e1; % reward high track centerline velocites
    q_r   = 0*1e5; % penalization when vehicle is outside track
    
    % state and input penalization
    q_br  = 0*1e2; % penalization of breaking
    q_acc = 0*1e2; % reward for acceleration
    q_v   = 0*1e1; % reward high absolute velocities

    % label inputs and outputs
    I_x = mu_x(1);          % x position in global coordinates
    I_y = mu_x(2);          % y position in global coordinates
    V_vx = mu_x(4);         % x velocity in vehicle coordinates
    psi  = mu_x(3);
    track_dist = mu_x(8);   % track velocity
    T = u(2);               % torque gain (1=max.acc, -1=max.braking)
    track_vel = u(4);       % track velocity
    
    % get information (x,y,track radius and track orientation) of the point 
    % in the track that corresponds to a traveled distance of 'dist' meters.
    [pos_c, psi_c, R_c] = track.getTrackInfo(track_dist);
    
    % ---------------------------------------------------------------------
    % cost of contour and lag error
    % ---------------------------------------------------------------------
    % rotation to a frame with x-axis tangencial to the track (T frame)
    A_TI = [cos(psi_c) -sin(psi_c);    
            sin(psi_c)  cos(psi_c)];
    % error in the inertial coordinates
    I_error = pos_c -[I_x;I_y];       
    % error in the T frame [lag_error; contouring_error]
    T_error = A_TI * I_error;          
    
    cost_contour = q_c*T_error(2)^2;
    cost_lag     = q_l*T_error(1)^2;
    
    % ---------------------------------------------------------------------
    % cost for orientation error (vehicle aligned with track orientation)
    % ---------------------------------------------------------------------
    cost_orientation = q_o*(psi_c-psi)^2;
    
    % ---------------------------------------------------------------------
    % cost for being outside track
    % ---------------------------------------------------------------------
    % is the vehicle outside the track?
    isOusideTrack = abs(T_error(2)) > R_c;
    if isOusideTrack
        % warning('Vehicle is outside the track!!!');
    end
    cost_outside = isOusideTrack*q_r*(norm(I_error)-R_c)^2;
    
    % ---------------------------------------------------------------------
    % reward high velocities only if inside track
    % ---------------------------------------------------------------------
    cost_vel = (~isOusideTrack)*-q_v*V_vx;
    
    % ---------------------------------------------------------------------
    % reward track velocities only if inside track
    % ---------------------------------------------------------------------
    cost_dist = (~isOusideTrack)*-q_d*track_vel;
    
    % ---------------------------------------------------------------------
    % reward acceleration and penalize braking
    % ---------------------------------------------------------------------
    cost_inputs = (T<0)*q_br*(T)^2 - (T>0)*q_acc*T^2;
    
    % ---------------------------------------------------------------------
    % Calculate final cost
    % ---------------------------------------------------------------------
    cost = cost_contour + cost_lag + cost_orientation + cost_dist + cost_outside + cost_inputs + cost_vel;
end


