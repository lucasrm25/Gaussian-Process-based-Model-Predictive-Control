%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot_racetrack
%
% plots the racetrack and your result
%
% files requested: racetrack.mat
%
% variables requested: Y (built by simulate_controlled_singletrack.m)
%
% The model is adopted from the scripts of M. Gerdts and D. Schramm,
% respectively.
%
% This file is for use within the "Project Competition" of the "Concepts of
% Automatic Control" course at the University of Stuttgart, held by F.
% Allgoewer.
%
% written by J. M. Montenbruck, Dec. 2013 
% mailto:jan-maximilian.montenbruck@ist.uni-stuttgart.de

load('racetrack.mat','t_r'); % load right  boundary from *.mat file
load('racetrack.mat','t_l'); % load left boundary from *.mat file
figure('Name','racetrack','NumberTitle','off','Toolbar','figure','MenuBar','none','OuterPosition',[0 -500 460 1100]) % creates window for plot
hold on % allow for multiple plot commands within one figure
axis equal % eqal axis scaling
axis([-50 70 -50 450]) % plot height and width
plot(t_r(:,1),t_r(:,2)) % plot right racetrack boundary
plot(t_l(:,1),t_l(:,2)) % plot left racetrack boundary
text(1,0,'\leftarrow finish/start','HorizontalAlignment','left') % finish/start annotation
plot(Y(:,1),Y(:,2),'r') % plot the x and y coordinates resulting fromy your controller
xlabel('x') % label x axis
ylabel('y') % label y axies
box % make a box around the plot