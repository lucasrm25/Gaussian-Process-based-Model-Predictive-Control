classdef RaceTrack
%------------------------------------------------------------------
%   object constructor
%------------------------------------------------------------------
    
    properties
    end
    
    properties(SetAccess=private)
        track_l         % <2,L>: left track coordinates
        track_r         % <2,L>: right track coordinates
        dist
    end
    
    methods
        function obj = RaceTrack(map)
        %------------------------------------------------------------------
        %   object constructor
        %------------------------------------------------------------------
            
        end
        
        
        function [closestPoint, latDist, longDist] = getClosestPointAhead(obj,x,distahead)
        %------------------------------------------------------------------
        %   
        %------------------------------------------------------------------
        end
        
        
        function generateTrackPoints(trackdata,x0, th0,w)
        %------------------------------------------------------------------
        %   
        %------------------------------------------------------------------
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
            obj.track_l = track_l;
            obj.track_r = track_r;
        end
    end
end

