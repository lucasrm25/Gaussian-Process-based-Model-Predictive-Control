classdef SingleTrackAnimation < handle
    
    properties
        racetrack @RaceTrack
        estPred  = NaN(4,11)       % <2,N>: approximated state predictions (x,y,vx,vy)
        ref      = NaN(2,11)       % <2,N>: references for each prediction
        truePred = NaN(4,11)       % <2,N>: true state predictions (x,y,vx,vy)
        
        figh
        ltrackh
        rtrackh
        estPredh
        refh
    end
    
    methods
        function obj = SingleTrackAnimation(racetrack, horizon)
            obj.racetrack = racetrack;
            obj.estPred  = NaN(4,horizon+1);
            obj.truePred = NaN(4,horizon+1);
        end
        
        function initGraphics(obj)
            obj.figh = figure('Color','w');
            hold on;
            grid on;
            axis equal;
            
            % -------------------------------------------------------------
            %   plot track asphalt
            % -------------------------------------------------------------
            n = length(obj.racetrack.track_l);
            v = [obj.racetrack.track_l(:,1:n)'; obj.racetrack.track_r(:,1:n)'];   % <n,2>
            f = [1:n-1 ; 2:n; n+2:2*n; n+1:2*n-1]';
            patch('Faces',f,'Vertices',v,'FaceColor',[0.95 0.95 0.95], 'LineStyle', 'none')
            
            % -------------------------------------------------------------
            %   plot track borders
            % -------------------------------------------------------------
            obj.ltrackh = plot(obj.racetrack.track_l(1,:),obj.racetrack.track_l(2,:));
            obj.rtrackh = plot(obj.racetrack.track_r(1,:),obj.racetrack.track_r(2,:));
            
            % -------------------------------------------------------------
            %   setup state predictions
            % -------------------------------------------------------------
            
            obj.refh = plot(obj.ref(1,:),obj.ref(2,:),...
                        '-','LineWidth',1.5, 'Marker','o',...
                        'MarkerFaceColor','b','MarkerEdgeColor','k', 'Color','b',...
                        'MarkerSize',5,...
                        'DisplayName','Projected optimal trajectory');
            
            obj.estPredh = patch(obj.estPred(1,:),obj.estPred(2,:),obj.estPred(3,:),...
                        'EdgeColor','interp','Marker','o',...
                        'MarkerSize',5,...
                        'MarkerFaceColor','flat',...
                        'DisplayName','Optimal trajectory');
            
%             obj.estPredh = plot(obj.estPred(1,:),obj.estPred(2,:),...
%                         '-','LineWidth',1, 'Marker','o',...
%                         'MarkerFaceColor','red','MarkerEdgeColor','k', 'Color','k',...
%                         'MarkerSize',5,...
%                         'DisplayName','Optimal trajectory');
  
            legend([obj.estPredh,obj.refh],'Location','southwest');
            c = colorbar;
            c.Label.String = 'Vehicle predicted velocity [m/s]';
            drawnow
        end
        
        function updateGraphics(obj)
            % update predicted trajectory
            obj.estPredh.XData = [obj.estPred(1,:) 0];
            obj.estPredh.YData = [obj.estPred(2,:) NaN];
            vel = vecnorm(obj.estPred(3:4,:));
            obj.estPredh.CData = [vel min(vel)];
            % update projected reference
            obj.refh.XData = obj.ref(1,:);
            obj.refh.YData = obj.ref(2,:);
            drawnow
        end
    end
end

