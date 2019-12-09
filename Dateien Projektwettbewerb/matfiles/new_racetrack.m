%% racetrack
%
% builds the racetrack and saves it as racetrack.mat
%
% files built: racetrack.mat 
%
% circle segment
% y = sqrt(r^2-(x-x0)^2)+ y0

% segment 1: straight, length 50
t_1_r=[zeros(50,1) (linspace(0,50,50))']; % right racetrack boundary, segment 1
t_1_l=[-5*ones(50,1) (linspace(0,50,50))']; % left racetrack boundary, segment 1

% segment 2: quater curve, inner radius 10
t_2_r=[(linspace(0,10,100))'  sqrt(10^2-((linspace(0,10,100))'-10).^2)+50]; % right racetrack boundary, segment 2
t_2_l=[(linspace(-5,10,100))' sqrt(15^2-((linspace(-5,10,100))'-10).^2)+50]; % left racetrack boundary, segment 2

% segment 3: straight, length 20
t_3_r=[(linspace(10,30,20))' 60*ones(20,1)]; % right racetrack boundary, segment 3
t_3_l=[(linspace(10,30,20))' 65*ones(20,1)]; % left racetrack boundary, segment 3
% segment 4: quater curve, inner radius 30
t_4_r=[sqrt(30^2-((linspace(60,30,80))'-30).^2)+40 (linspace(60,30,80))'  ]; % right racetrack boundary, segment 4
t_4_l=[sqrt(35^2-((linspace(65,30,80))'-30).^2)+40 (linspace(65,30,80))' ]; % left racetrack boundary, segment 4

% segment 5: straight, length 5
t_5_r=[70*ones(10,1) (linspace(30,25,10))']; % right racetrack boundary, segment 5
t_5_l=[75*ones(10,1) (linspace(30,25,10))']; % left racetrack boundary, segment 5

% segment 6: quarter curve, inner radius 2.5
t_6_r=[(linspace(70,77.5,20))' -sqrt(7.5^2-((linspace(70,77.7,20))'-77.5).^2)+25]; % right racetrack boundary, segment 6
t_6_l=[(linspace(75,77.7,20))' -sqrt(2.5^2-((linspace(75,77.7,20))'-77.5).^2)+25]; % left racetrack boundary, segment 6

% segment 7: quarter curve, inner radius 2.5

t_7_r=[sqrt(2.5^2-((linspace(17.5,15,20))'-15).^2)+77.5 (linspace(17.5,15,20))']; % right racetrack boundary, segment 7
t_7_l=[sqrt(7.5^2-((linspace(22.5,15,20))'-15).^2)+77.5 (linspace(22.5,15,20))']; % left racetrack boundary, segment 7

% segment 8: straight, length 10
t_8_r=[80*ones(20,1) (linspace(15,5,20))']; % right racetrack boundary, segment 8
t_8_l=[85*ones(20,1) (linspace(15,5,20))']; % left racetrack boundary, segment 8

% segment 9: halv curve, inner radius 5
t_9_r=[(linspace(80,70,50))' -sqrt(5^2-((linspace(80,70,50))'-75).^2)+5]; % right racetrack boundary, segment 9
t_9_l=[(linspace(85,65,50))' -sqrt(10^2-((linspace(85,65,50))'-75).^2)+5]; % left racetrack boundary, segment 9

% segment 10: halv curve, inner radius 5
t_10_r=[(linspace(70,50,50))' sqrt(10^2-((linspace(70,50,50))'-60).^2)+5]; % right racetrack boundary, segment 10
t_10_l=[(linspace(65,55,50))' sqrt(5^2-((linspace(65,55,50))'-60).^2)+5]; % left racetrack boundary, first segment 10


% segment 11: straight, length 5
t_11_r=[50*ones(20,1) (linspace(5,0,20))']; % right racetrack boundary, segment 11
t_11_l=[55*ones(20,1) (linspace(5,0,20))']; % left racetrack boundary, segment 11

% segment 12: halv curve, inner radius 25
t_12_r=[(linspace(50,0,100))' -sqrt(25^2-((linspace(50,0,100))'-25).^2)+0]; % right racetrack boundary, segment 12
t_12_l=[(linspace(55,-5,100))' -sqrt(30^2-((linspace(55,-5,100))'-25).^2)+0]; % left racetrack boundary, segment 12
%% 
t_r=[t_1_r ; t_2_r ; t_3_r ; t_4_r ; t_5_r ; t_6_r ; t_7_r ; t_8_r ; t_9_r ; t_10_r ; t_11_r ; t_12_r ];%; t_13_r ; t_14_r]; % stack of right racetrack boundaries
t_l=[t_1_l ; t_2_l ; t_3_l ; t_4_l ; t_5_l ; t_6_l ; t_7_l ; t_8_l ; t_9_l ; t_10_l ; t_11_l ; t_12_l ];% ; t_13_l ; t_14_l]; % stack of left racetrack boundaries
save('racetrack.mat', 't_r', 't_l'); % save racetrack as racetrack.mat