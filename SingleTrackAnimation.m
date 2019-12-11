classdef SingleTrackAnimation < handle
    
    properties
        racetrack @RaceTrack
        truePred = NaN(2,11)       % <2,N>: true state predictions
        estPred  = NaN(2,11)       % <2,N>: approximated state predictions
        ref      = NaN(2,11)       % <2,N>: references for each prediction
        
        figh
        ltrackh
        rtrackh
        estPredh
        refh
    end
    
    methods
        function obj = SingleTrackAnimation(racetrack, horizon)
            obj.racetrack = racetrack;
            obj.truePred = NaN(2,horizon+1);
            obj.estPred  = NaN(2,horizon+1);
        end
        
        function initGraphics(obj)
            obj.figh = figure;
            hold on;
            grid on;
            axis equal;
            
            obj.ltrackh = plot(obj.racetrack.track_l(1,:),obj.racetrack.track_l(2,:));
            obj.rtrackh = plot(obj.racetrack.track_r(1,:),obj.racetrack.track_r(2,:));
            
            obj.estPredh = plot(obj.estPred(1,:),obj.estPred(2,:),'-*','LineWidth',1, 'DisplayName','Estimated Prediction');
            obj.refh = plot(obj.ref(1,:),obj.ref(2,:),'-*','LineWidth',1, 'DisplayName','Optimized target trajectory');
            drawnow
        end
        
        function updateGraphics(obj)
            obj.estPredh.XData = obj.estPred(1,:);
            obj.estPredh.YData = obj.estPred(2,:);
            obj.refh.XData = obj.ref(1,:);
            obj.refh.YData = obj.ref(2,:);
            drawnow
        end
    end
end

