%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% racetrack
%
% builds the racetrack and saves it as racetrack.mat


% t_1_r=[zeros(250,1) (linspace(0,250,250))']; % right racetrack boundary, segment 1
% t_1_l=[-5*ones(250,1) (linspace(0,250,250))']; % left racetrack boundary, segment 1
% t_2_r=[(linspace(0,-20,100))' sqrt(100-((linspace(0,-20,100))'+10).^2)+250]; % right racetrack boundary, segment 2
% t_2_l=[(linspace(-5,-15,100))' sqrt(25-((linspace(-5,-15,100))'+10).^2)+250]; % left racetrack boundary, segment 2
% t_3_r=[(linspace(-20,-30,100))' -sqrt(25-((linspace(-20,-30,100))'+25).^2)+250]; % right racetrack boundary, segment 3
% t_3_l=[(linspace(-15,-35,100))' -sqrt(100-((linspace(-15,-35,100))'+25).^2)+250]; % left racetrack boundary, segment 3
% t_4_r=[-30*ones(150,1) (linspace(250,400,150))']; % right racetrack boundary, segment 4
% t_4_l=[-35*ones(150,1) (linspace(250,400,150))']; % left racetrack boundary, segment 4
% t_5_r=[(linspace(-30,20,200))' sqrt(625-((linspace(-30,20,200))'+5).^2)+400]; % right racetrack boundary, segment 5
% t_5_l=[(linspace(-35,25,200))' sqrt(900-((linspace(-35,25,200))'+5).^2)+400]; % left racetrack boundary, segment 5
% t_6_r=[20*ones(100,1) (linspace(400,300,100))']; % right racetrack boundary, segment 6
% t_6_l=[25*ones(95,1) (linspace(400,305,95))']; % left racetrack boundary, segment 6
% t_7_r=[(linspace(20,25,50))' -sqrt(25-((linspace(20,25,50))'-25).^2)+300]; % right racetrack boundary, segment 7
% t_7_l=[(linspace(25,30,50))' -sqrt(25-((linspace(25,30,50))'-30).^2)+305]; % left racetrack boundary, segment 7
% t_8_r=[(linspace(25,30,50))' sqrt(25-((linspace(25,30,50))'-25).^2)+290]; % right racetrack boundary, segment 8
% t_8_l=[(linspace(30,35,50))' sqrt(25-((linspace(30,35,50))'-30).^2)+295]; % left racetrack boundary, segment 8
% t_9_r=[30*ones(245,1) (linspace(290,45,245))']; % right racetrack boundary, segment 9
% t_9_l=[35*ones(250,1) (linspace(295,45,250))']; % left racetrack boundary, segment 9
% t_10_r=[(linspace(30,45,150))' -sqrt(225-((linspace(30,45,150))'-45).^2)+45]; % right racetrack boundary, segment 10
% t_10_l=[(linspace(35,45,150))' -sqrt(100-((linspace(35,45,150))'-45).^2)+45]; % left racetrack boundary, first segment 10
% t_11_r=[(linspace(45,55,150))' sqrt(100-((linspace(45,55,150))'-45).^2)+20]; % right racetrack boundary, segment 11
% t_11_l=[(linspace(45,60,150))' sqrt(225-((linspace(45,60,150))'-45).^2)+20]; % left racetrack boundary, segment 11
% t_12_r=[(linspace(55,20,250))' -sqrt(1225-((linspace(55,20,250))'-20).^2)+20]; % right racetrack boundary, segment 12
% t_12_l=[(linspace(60,20,250))' -sqrt(1600-((linspace(60,20,250))'-20).^2)+20]; % left racetrack boundary, segment 12
% t_13_r=[(linspace(20,15,5))' -15*ones(5,1)]; % right racetrack boundary, segment 13
% t_13_l=[(linspace(20,15,5))' -20*ones(5,1)]; % left racetrack boundary, segment 13
% t_14_r=[(linspace(15,0,150))' -sqrt(225-((linspace(15,0,150))'-15).^2)]; % right racetrack boundary, segment 14
% t_14_l=[(linspace(15,-5,150))' -sqrt(400-((linspace(15,-5,150))'-15).^2)]; % left racetrack boundary, segment 14
% t_r=[t_1_r ; t_2_r ; t_3_r ; t_4_r ; t_5_r ; t_6_r ; t_7_r ; t_8_r ; t_9_r ; t_10_r ; t_11_r ; t_12_r ; t_13_r ; t_14_r]; % stack of right racetrack boundaries
% t_l=[t_1_l ; t_2_l ; t_3_l ; t_4_l ; t_5_l ; t_6_l ; t_7_l ; t_8_l ; t_9_l ; t_10_l ; t_11_l ; t_12_l ; t_13_l ; t_14_l]; % stack of left racetrack boundaries


%% New track

track_nbr = 2;


if track_nbr == 1
    x0  = [0;0];
    th0 = 0;
    w = 3;
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
else
    x0  = [0;0];
    th0 = 0;
    w = 2;
    trackdata = {
     's',17;
     'c',{8,-90};
     's',5;
     'c',{2,90};
     'c',{2,-90};
     's',12;
     'c',{2.2,-180};
     's',1;
     'c',{2.2,180};
     's',1;
     'c',{2.2,-90};
     's',10;
     'c',{2.2,-140};
     's',12;
     'c',{8,140};
     's',6;
     'c',{8,90};
     's',3;
     'c',{2.5,130};
     's',10;
     'c',{2.2,-180};
     's',22;
     'c',{2.2,-130};
     's',17.8;      
     'c',{10,-90}; 
    };
end
    
    
[t_l,t_r] = calc_track(trackdata,x0,th0,w);


% fig = figure; hold on; grid on; axis equal;
% plot(t_l(1,:),t_l(2,:))
% plot(t_r(1,:),t_r(2,:));
% pause;
% close(fig);


%% Save track

assert( size(t_r,2)==2 && size(t_l,2)==2, 'track must be a Nx2 vector');
save('racetrack.mat', 't_r', 't_l'); % save racetrack as racetrack.mat



%% Calculate track patches

function [track_l,track_r] = calc_track(trackdata,x0, th0,w)
    track_l = [];
    track_r = [];
    ds = 0.5;
    dth = deg2rad(2);
    
    A_z = @(th) [cos(th) -sin(th);
                 sin(th)  cos(th)];
             
    A_IV = A_z(th0);
    r_IV = x0;
    
    for idx = 1:size(trackdata,1)
        if strcmp( trackdata{idx,1},'s')
            % distance
            dist = trackdata{idx,2};
            % calculate new track
            newtrack_l = r_IV + w*A_IV(:,2) + (ds:ds:dist).*A_IV(:,1);
            newtrack_r = r_IV - w*A_IV(:,2) + (ds:ds:dist).*A_IV(:,1);
            % update current position
            r_IV = r_IV + dist*A_IV(:,1);
            
        elseif strcmp( trackdata{idx,1},'c')
            % center curve radius
            rad = trackdata{idx,2}{1};
            % curvature
            ang = deg2rad( trackdata{idx,2}{2} );
            
            th = 0:dth:abs(ang);
            arc   = [cos(th); sin(th)];
            if ang > 0
                newtrack_l = arc*(rad-w);
                newtrack_r = arc*(rad+w);
                A = [ 0 1 
                     -1 0];
                newtrack_l = A*newtrack_l + [0;1]*rad;
                newtrack_r = A*newtrack_r + [0;1]*rad; 
            else
                newtrack_l = arc*(rad+w);
                newtrack_r = arc*(rad-w);
                A = [0 1 
                     1 0];
                newtrack_l = A*newtrack_l + [0;-1]*rad;
                newtrack_r = A*newtrack_r + [0;-1]*rad; 
            end

            newtrack_l = A_IV * newtrack_l + r_IV;
            newtrack_r = A_IV * newtrack_r + r_IV;
            
            A_IV = A_z(ang) * A_IV;
            r_IV = 0.5* (newtrack_l(:,end) + newtrack_r(:,end));
        else
            error('Not implemented');
        end
        
        track_l = [track_l newtrack_l];
        track_r = [track_r newtrack_r];
    end
    track_l = track_l';
    track_r = track_r';
end

