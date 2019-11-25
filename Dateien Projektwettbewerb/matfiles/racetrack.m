%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% racetrack
%
% builds the racetrack and saves it as racetrack.mat
%
% files built: racetrack.mat 
%
% This file is for use within the "Project Competition" of the "Concepts of
% Automatic Control" course at the University of Stuttgart, held by F.
% Allgoewer.
%
% written by J. M. Montenbruck, Dec. 2013 
% mailto:jan-maximilian.montenbruck@ist.uni-stuttgart.de

t_1_r=[zeros(250,1) (linspace(0,250,250))']; % right racetrack boundary, segment 1
t_1_l=[-5*ones(250,1) (linspace(0,250,250))']; % left racetrack boundary, segment 1
t_2_r=[(linspace(0,-20,100))' sqrt(100-((linspace(0,-20,100))'+10).^2)+250]; % right racetrack boundary, segment 2
t_2_l=[(linspace(-5,-15,100))' sqrt(25-((linspace(-5,-15,100))'+10).^2)+250]; % left racetrack boundary, segment 2
t_3_r=[(linspace(-20,-30,100))' -sqrt(25-((linspace(-20,-30,100))'+25).^2)+250]; % right racetrack boundary, segment 3
t_3_l=[(linspace(-15,-35,100))' -sqrt(100-((linspace(-15,-35,100))'+25).^2)+250]; % left racetrack boundary, segment 3
t_4_r=[-30*ones(150,1) (linspace(250,400,150))']; % right racetrack boundary, segment 4
t_4_l=[-35*ones(150,1) (linspace(250,400,150))']; % left racetrack boundary, segment 4
t_5_r=[(linspace(-30,20,200))' sqrt(625-((linspace(-30,20,200))'+5).^2)+400]; % right racetrack boundary, segment 5
t_5_l=[(linspace(-35,25,200))' sqrt(900-((linspace(-35,25,200))'+5).^2)+400]; % left racetrack boundary, segment 5
t_6_r=[20*ones(100,1) (linspace(400,300,100))']; % right racetrack boundary, segment 6
t_6_l=[25*ones(95,1) (linspace(400,305,95))']; % left racetrack boundary, segment 6
t_7_r=[(linspace(20,25,50))' -sqrt(25-((linspace(20,25,50))'-25).^2)+300]; % right racetrack boundary, segment 7
t_7_l=[(linspace(25,30,50))' -sqrt(25-((linspace(25,30,50))'-30).^2)+305]; % left racetrack boundary, segment 7
t_8_r=[(linspace(25,30,50))' sqrt(25-((linspace(25,30,50))'-25).^2)+290]; % right racetrack boundary, segment 8
t_8_l=[(linspace(30,35,50))' sqrt(25-((linspace(30,35,50))'-30).^2)+295]; % left racetrack boundary, segment 8
t_9_r=[30*ones(245,1) (linspace(290,45,245))']; % right racetrack boundary, segment 9
t_9_l=[35*ones(250,1) (linspace(295,45,250))']; % left racetrack boundary, segment 9
t_10_r=[(linspace(30,45,150))' -sqrt(225-((linspace(30,45,150))'-45).^2)+45]; % right racetrack boundary, segment 10
t_10_l=[(linspace(35,45,150))' -sqrt(100-((linspace(35,45,150))'-45).^2)+45]; % left racetrack boundary, first segment 10
t_11_r=[(linspace(45,55,150))' sqrt(100-((linspace(45,55,150))'-45).^2)+20]; % right racetrack boundary, segment 11
t_11_l=[(linspace(45,60,150))' sqrt(225-((linspace(45,60,150))'-45).^2)+20]; % left racetrack boundary, segment 11
t_12_r=[(linspace(55,20,250))' -sqrt(1225-((linspace(55,20,250))'-20).^2)+20]; % right racetrack boundary, segment 12
t_12_l=[(linspace(60,20,250))' -sqrt(1600-((linspace(60,20,250))'-20).^2)+20]; % left racetrack boundary, segment 12
t_13_r=[(linspace(20,15,5))' -15*ones(5,1)]; % right racetrack boundary, segment 13
t_13_l=[(linspace(20,15,5))' -20*ones(5,1)]; % left racetrack boundary, segment 13
t_14_r=[(linspace(15,0,150))' -sqrt(225-((linspace(15,0,150))'-15).^2)]; % right racetrack boundary, segment 14
t_14_l=[(linspace(15,-5,150))' -sqrt(400-((linspace(15,-5,150))'-15).^2)]; % left racetrack boundary, segment 14
t_r=[t_1_r ; t_2_r ; t_3_r ; t_4_r ; t_5_r ; t_6_r ; t_7_r ; t_8_r ; t_9_r ; t_10_r ; t_11_r ; t_12_r ; t_13_r ; t_14_r]; % stack of right racetrack boundaries
t_l=[t_1_l ; t_2_l ; t_3_l ; t_4_l ; t_5_l ; t_6_l ; t_7_l ; t_8_l ; t_9_l ; t_10_l ; t_11_l ; t_12_l ; t_13_l ; t_14_l]; % stack of left racetrack boundaries
save('racetrack.mat', 't_r', 't_l'); % save racetrack as racetrack.mat