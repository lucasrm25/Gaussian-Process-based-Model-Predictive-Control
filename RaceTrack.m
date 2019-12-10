classdef RaceTrack
%------------------------------------------------------------------
%   The Race Track constitutes of a left lane, right lane and center line
%   coordinates, which are parametrized by the traveled distance. 
%   
%   This also means that the center line position, pos_c = (x_c,y_c), the
%   track radius and the track orientation can be obtained as a function of
%   the traveled distance:
%
%   x_c(dist), y_c(dist), psi_c(dist), R_c(dist) -> see method getTrackInfo
%   
%   
%   To create the track, this Class requires a series of instructions. 
%    Ex:
%
%         x0  = [0;0];              % initial track position
%         th0 = 0;                  % initial track orientation
%         w = 6;                    % track width
%         trackdata = {
%              's',14;              % go straight ahead 14 meters
%              'c',{15,-90};        % make a turn of 90deg and radius 15 meters
%              's',5;               % go straight ahead 5 meters
%              'c',{4,90};          % ...
%              'c',{4,-90};
%              's',5;
%              'c',{3.5,-180};
%              'c',{3.5,180};
%              'c',{3.5,-90};
%              's',2;
%              'c',{3.5,-120};
%              's',10;
%              'c',{10,120};
%              's',10;
%              'c',{5,90};
%              's',5;
%              'c',{5,150};
%              's',5;
%              'c',{3.2,-180};
%              's',12;
%              'c',{10,-150};
%              's',12.3;      
%              'c',{12,-90}; 
%         };
%         track = RaceTrack(trackdata, x0, th0, w);
%------------------------------------------------------------------
    
    properties
    end
    
    properties(SetAccess=private)
        w
        track_l         % <2,L>: left track coordinates
        track_r         % <2,L>: right track coordinates
        
        track_c
        psi_c
        R_c
        
        dist            % traveled distance from the starting point
    end
    
    methods
        function obj = RaceTrack(trackdata, x0, th0, w)
        %------------------------------------------------------------------
        %   object constructor
        %   args:
        %       x0: <2,1> initial vehicle position
        %       tho: <1> initial vehicle orientation
        %       w: <1> track width
        %       trackdata: sequence of track instructions
        %------------------------------------------------------------------
            [obj.track_l, obj.track_r, obj.psi_c] = obj.generateTrackPoints(trackdata, x0, th0, w);
            obj.track_c = (obj.track_r + obj.track_l)/2;
            
            % parametrize track by traveled distance
            distFromCenter = sqrt(obj.track_c(1,:).^2 + obj.track_c(2,:).^2);
            obj.dist = cumsum( abs(conv(distFromCenter,[1+1e-6,-1],'valid')) );
            
            obj.w = w;
        end
        
        
        function [pos_c, psi_c, R_c] = getTrackInfo(obj,dist)
        %------------------------------------------------------------------
        %   Given a distance 'dist' returns the track centerline carthesian 
        %   position (Xt,Yt), the track orientation PSIt, and the track
        %   radius Rt
        %------------------------------------------------------------------
            dist = mod(dist,max(obj.dist));
            idx_dist = interp1(obj.dist,1:numel(obj.dist),dist,'nearest','extrap');
            
            x_c   = obj.track_c(1,idx_dist);
            y_c   = obj.track_c(2,idx_dist);
            pos_c = [x_c;y_c];
            psi_c = obj.psi_c(1,idx_dist);
            R_c   = obj.w/2;
        end
    end
       
    methods(Static)
        function [track_l,track_r,psi_c] = generateTrackPoints(trackdata, x0, th0, w)
        %------------------------------------------------------------------
        %   
        %------------------------------------------------------------------
            track_l = [];
            track_r = [];
            psi_c = [];
            
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
                    newtrack_l = r_IV + w/2*A_IV(:,2) + (ds:ds:dist).*A_IV(:,1);
                    newtrack_r = r_IV - w/2*A_IV(:,2) + (ds:ds:dist).*A_IV(:,1);
                    newpsi_c   = th0 * ones(1,size(newtrack_l,2));
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
                        newtrack_l = arc*(rad-w/2);
                        newtrack_r = arc*(rad+w/2);
                        A = [ 0 1 
                             -1 0];
                        newtrack_l = A*newtrack_l + [0;1]*rad;
                        newtrack_r = A*newtrack_r + [0;1]*rad; 
                    else
                        newtrack_l = arc*(rad+w/2);
                        newtrack_r = arc*(rad-w/2);
                        A = [0 1 
                             1 0];
                        newtrack_l = A*newtrack_l + [0;-1]*rad;
                        newtrack_r = A*newtrack_r + [0;-1]*rad; 
                    end

                    newtrack_l = A_IV * newtrack_l + r_IV;
                    newtrack_r = A_IV * newtrack_r + r_IV;
                    newpsi_c   = th0 + th;

                    A_IV = A_z(ang) * A_IV;
                    r_IV = 0.5* (newtrack_l(:,end) + newtrack_r(:,end));
                else
                    error('Not implemented');
                end

                track_l = [track_l newtrack_l];
                track_r = [track_r newtrack_r];
                psi_c   = [psi_c newpsi_c];
            end
        end
    end
end

