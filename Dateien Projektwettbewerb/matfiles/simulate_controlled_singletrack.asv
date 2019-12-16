function simulate_controlled_singletrack(t_f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function simulate_controlled_singletrack(t_f)
%
% integrates the controlled single-track model until time t_f
%
% input: t_f (simulation time)
%
% files requested: racetrack.m ; singletrack.m ; ode1.m ; plot_racetrack.m
%
% plots built: racetrack
%
% This file is for use within the "Project Competition" of the "Concepts of
% Automatic Control" course at the University of Stuttgart, held by F.
% Allgoewer.
%
% written by J. M. Montenbruck, Dec. 2013 
% mailto:jan-maximilian.montenbruck@ist.uni-stuttgart.de

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_racetrack % builds the racetrack and saves it as racetrack.mat


dt = 0.1;



% -------------------------------------------------------------------------

  % load a predefined track called track01:
      [trackdata, x0, th0, w] = RaceTrack2.loadTrack_01();

  % init object:
      global trackr;
      trackr = RaceTrack2(trackdata, x0, th0, w);
      




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTEGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X_0=[-2.5;0;0;0;pi/2;0;0;0;0;0]; % initial value for integration
%X_0=[72.54;31.38;7.5;0;3/2*pi;0;0;0;0;0]; 
X_0= [x0(1);x0(2);0;0;0;0;0;0;0;0];
Y=ode1(@singletrack,0:dt:t_f,X_0); % integrate with step zise 0.001

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EVALUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_new_racetrack % plots the racetrack and your result
end