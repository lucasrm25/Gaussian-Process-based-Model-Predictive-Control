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
t_2_r=[(linspace(0,10,20))'  sqrt(10^2-((linspace(0,10,20))'-10).^2)+50]; % right racetrack boundary, segment 2
t_2_l=[(linspace(-5,10,20))' sqrt(15^2-((linspace(-5,10,20))'-10).^2)+50]; % left racetrack boundary, segment 2

% segment 3: straight, length 20
t_3_r=[(linspace(10,30,30))' 60*ones(30,1)]; % right racetrack boundary, segment 3
t_3_l=[(linspace(10,30,30))' 65*ones(30,1)]; % left racetrack boundary, segment 3

% segment 4: quater curve, inner radius 30
t_4_r=[sqrt(30^2-((linspace(60,30,80))'-30).^2)+40 (linspace(60,30,80))'  ]; % right racetrack boundary, segment 4
t_4_l=[sqrt(35^2-((linspace(65,30,80))'-30).^2)+40 (linspace(65,30,80))' ]; % left racetrack boundary, segment 4

% segment 5: straight, length 5
t_5_r=[70*ones(5,1) (linspace(30,25,5))']; % right racetrack boundary, segment 5
t_5_l=[75*ones(5,1) (linspace(30,25,5))']; % left racetrack boundary, segment 5

% segment 6: quarter curve, inner radius 2.5
t_6_r=[(linspace(70,77.5,5))' -sqrt(7.5^2-((linspace(70,77.5,5))'-77.5).^2)+25]; % right racetrack boundary, segment 6
t_6_l=[(linspace(75,77.5,5))' -sqrt(2.5^2-((linspace(75,77.5,5))'-77.5).^2)+25]; % left racetrack boundary, segment 6

% segment 11: straight, length 5
t_11_r=[(linspace(77.5,82.5,15))' 17.5*ones(15,1) ]; % right racetrack boundary, segment 11
t_11_l=[(linspace(77.5,82.5,15))' 22.5*ones(15,1) ]; % left racetrack boundary, segment 11

% segment 7: quarter curve, inner radius 2.5
% t_7_r=[sqrt(2.5^2-((linspace(17.5,15,5))'-15).^2)+82.5 (linspace(17.5,15,5))']; % right racetrack boundary, segment 7
% t_7_l=[sqrt(7.5^2-((linspace(22.5,15,5))'-15).^2)+82.5 (linspace(22.5,15,5))']; % left racetrack boundary, segment 7

% segment 8: straight, length 10
% t_8_r=[85*ones(15,1) (linspace(15,0,15))']; % right racetrack boundary, segment 8
% t_8_l=[90*ones(15,1) (linspace(15,0,15))']; % left racetrack boundary, segment 8

t_8_r=[sqrt(12.5^2-((linspace(17.5,-7.5,40))'-5).^2)+82.5 (linspace(17.5,-7.5,40))' ]; % right racetrack boundary, segment 10
t_8_l=[sqrt(17.5^2-((linspace(22.5,-12.5,40))'-5).^2)+82.5 (linspace(22.5,-12.5,40))' ]; % left racetrack boundary, first segment 10

t_91_r=[(linspace(82.5,77.5,5))' -7.5*ones(5,1)]; % right racetrack boundary, segment 9
t_91_l=[(linspace(82.5,77.5,5))' -12.5*ones(5,1)]; % left racetrack boundary, segment 9

% segment 9: halv curve, inner radius 5
t_9_r=[(linspace(77.5,70,30))' -sqrt(7.5^2-((linspace(77.5,70,30))'-77.5).^2)+0]; % right racetrack boundary, segment 9
t_9_l=[(linspace(77.5,65,30))' -sqrt(12.5^2-((linspace(77.5,65,30))'-77.5).^2)+0]; % left racetrack boundary, segment 9

% segment 10: halv curve, inner radius 5
t_10_r=[(linspace(70,50,20))' sqrt(10^2-((linspace(70,50,20))'-60).^2)+0]; % right racetrack boundary, segment 10
t_10_l=[(linspace(65,55,20))' sqrt(5^2-((linspace(65,55,20))'-60).^2)+0]; % left racetrack boundary, first segment 10



% segment 12: halv curve, inner radius 25
t_12_r=[(linspace(50,0,150))' -sqrt(25^2-((linspace(50,0,150))'-25).^2)+0]; % right racetrack boundary, segment 12
t_12_l=[(linspace(55,-5,150))' -sqrt(30^2-((linspace(55,-5,150))'-25).^2)+0]; % left racetrack boundary, segment 12
%% 
t_r=[t_1_r ; t_2_r ; t_3_r ; t_4_r ; t_5_r ; t_6_r ;t_11_r ; t_8_r ;t_91_r ; t_9_r ; t_10_r  ; t_12_r ];%; t_13_r ; t_14_r]; % stack of right racetrack boundaries
t_l=[t_1_l ; t_2_l ; t_3_l ; t_4_l ; t_5_l ; t_6_l ;t_11_l ; t_8_l ; t_91_l; t_9_l ; t_10_l  ; t_12_l ];% ; t_13_l ; t_14_l]; % stack of left racetrack boundaries

t_c = 0.5 * (t_r + t_l);


distmax = 0.4;
new_t_c = [];

for k=1:length(t_c)-1
    n = norm( t_c(k,:) - t_c(k+1,:) ) / distmax;
    i = 1;
    for alpha = (1:n-1e-3)/n
        newpoints(i,:) = alpha * t_c(k+1,:) + (1-alpha) * t_c(k,:);
        i = i+1;
    end
%     new_t_c = [new_t_c(1:k,:); newpoints; new_t_c(k+1:end,:) ];
    new_t_c = [new_t_c; t_c(k,:); newpoints];
    newpoints = [];
end

% figure
% plot(new_t_c(1:end,1),new_t_c(1:end,2),'.')
% 
% figure
% plot(t_c(1:450,1),t_c(1:450,2),'.')

save('racetrack.mat', 't_r', 't_l', 'new_t_c'); % save racetrack as racetrack.mat