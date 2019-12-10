classdef SingleTrackAnimation
    
    properties
        racetrack       @RaceTrack
        TruePred        % <2,N>: true state predictions
        AppPred         % <2,N>: approximated state predictions
        ref             % <2,N>: references for each prediction
        
        figh
        ltrackh
        rtrackh
    end
    
    methods
        function obj = SingleTrackAnimation(racetrack)
            obj.racetrack = racetrack;
        end
        
        function initGraphics(obj)
            obj.figh = figure;
            hold on;
            grid on;
            axis equal;
            
            obj.ltrackh = plot(obj.racetrack.track_l(1,:),obj.racetrack.track_l(2,:));
            obj.rtrackh = plot(obj.racetrack.track_r(1,:),obj.racetrack.track_r(2,:));
        end
        
        function updateGraphics(obj)
        end
    end
end

