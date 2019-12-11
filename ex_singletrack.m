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


% -------------------------------------------------------------------------
% Nonlinear Model Predictive Controller
% -------------------------------------------------------------------------

% define cost function
n  = estModel.n;
m  = estModel.m;
ne = 0;
N = 10; % prediction horizon
Q = diag([1000 1000 100]);
Qf= Q;
R = diag([0 0 0.1 0 -1]);
Ck = [eye(3), zeros(3,7)];

% define cost functions
fo   = @(t,mu_x,var_x,u,e,r) costFunction(mu_x, var_x, u, track);            % e = track distance
fend = @(t,mu_x,var_x,e,r)   costFunction(mu_x, var_x, zeros(m,1), track);   % end cost function

% define dynamics
f  = @(mu_x,var_x,u) estModel.xkp1(mu_x, var_x, u, dt);
% define additional constraints
h  = @(x,u,e) [];
g  = @(x,u,e) [];
u_lb = [-deg2rad(10);  % delta_dot >= -10deg/s
         1;            % gear  >= 1
         0;            % break >= 0
         0;            % break distribution >= 0
         0;            % acceleration pedal >= 0
         0];           % track velocity >= 0
u_ub = [deg2rad(10);   % delta_dot <=  10 deg/s
        5;             % gear  <= 5
        5000;          % break <= 5000
        1;             % break distribution <= 1
        1;             % acceleration pedal <= 1
        10];           % track velocity <= 1

% Initialize NMPC object;
mpc = NMPC(f, h, g, u_lb, u_ub, n, m, ne, fo, fend, N, dt);
mpc.tol     = 1e-5;
mpc.maxiter = 30;

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
true_x0(end) = track.getTrackDistance(true_x0(1:2));

est_x0  = [10;0;0; 5;0;0; 0; 0];   % initial state prior
% get initial track traveled distance
est_x0(end) = track.getTrackDistance(est_x0(1:2));  

% define simulation time
out.t = 0:dt:tf;            % time vector
kmax = length(out.t)-1;     % steps to simulate

% initialize variables to store simulation results
out.x    = [true_x0 NaN(true_n,kmax)];
out.xhat = [est_x0  NaN(est_n,kmax)];
out.xnom = [est_x0  NaN(est_n,kmax)];
out.u    = NaN(est_m,kmax);
% out.r = zeros(nr,length(out.t)-1);


% deactivate GP evaluation in the prediction
d_GP.isActive = false;

% start animation
trackAnim = SingleTrackAnimation(track,mpc.N);
trackAnim.initGraphics()

% scopex = figure('Position',[-1006 86 957 808]);
% scopeu = figure('Position',[-1879 93 867 795]);

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
    sprintf('\nSteering angle: %d\nGear: %d\nF_b: %.1f\nBrake dist: %.1f\nAcc pedal: %.1f\n',rad2deg(out.u(1,k)),out.u(2,k),out.u(3,k),out.u(4,k),out.u(5,k))
    
    % ---------------------------------------------------------------------
    % plot
    % ---------------------------------------------------------------------
    estPred = mpc.predictStateSequence(out.xhat(:,k), zeros(estModel.n), u_opt);
    trackAnim.estPred  = estPred(1:2,:); % we only need X,Y
    trackAnim.ref = track.getTrackInfo(estPred(end,:));
    trackAnim.updateGraphics();
    % plotScope(scopex,scopeu,out.xhat,out.u);
    
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


%% Evaluate results
close all;



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
    for i=1:8
        subplot(5,2,i);
        if angles(i)
            plot(rad2deg(x(i,:)),'DisplayName',names{i});
        else
            plot(x(i,:),'DisplayName',names{i});
        end
        legend('Location', 'west');
    end
    
    figure(figu);
    names = {'deltadot','G','F_b','Breaking dist.','Acc pedal'};
    angles = [1 0 0 0 0];
    for i=1:5
        subplot(5,1,i);
        if angles(i)
            plot(rad2deg(u(i,:)),'DisplayName',names{i});
        else
            plot(u(i,:),'DisplayName',names{i});
        end
        legend('Location', 'west');
    end
end

function cost = costFunction(mu_x, var_x, u, track)
    q_l   = 1e2; % penalization of lag error
    q_c   = 1e5; % penalization of contouring error
    q_r   = 1e1; % penalization when vehicle is outside track
    q_br  = 1e0; % penalization of breaking
    q_acc = 1e2; % reward for acceleration
    q_v   = 1e1; % reward high absolute velocities
    q_d   = 1e1; % reward high track centerline velocites

    % get information (x,y,track radius and track orientation) of the point 
    % in the track that corresponds to a traveled distance of 'dist' meters.
    track_dist = mu_x(8);
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
    isOusideTrack = abs(T_error(2)) > R_c;
    if isOusideTrack
        % warning('Vehicle is outside the track!!!');
    end
    cost_outside = isOusideTrack*q_r*(norm(I_error)-R_c)^2;
    
    % ---------------------------------------------------------------------
    % reward high velocities only if inside track
    % ---------------------------------------------------------------------
    cost_vel = (~isOusideTrack)*-q_v*mu_x(4);
    
    % ---------------------------------------------------------------------
    % reward track velocities only if inside track
    % ---------------------------------------------------------------------
    track_vel = u(6);
    cost_dist = (~isOusideTrack)*-q_d*track_vel;
    
    
    % ---------------------------------------------------------------------
    % cost for high inputs - penalize high inputs
    % ---------------------------------------------------------------------
    brake_force = u(3);
    acc_pedal   = u(5);
    cost_inputs = q_br*(brake_force/5000)^2 - q_acc*acc_pedal^2;
    
    % ---------------------------------------------------------------------
    % Calculate final cost
    % ---------------------------------------------------------------------
    cost = cost_contour + cost_lag + cost_outside + cost_inputs + cost_vel + cost_dist;
end


