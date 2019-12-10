classdef SingleTrackAnimation
    
    properties
        track_l         % <2,L>: left track coordinates
        track_r         % <2,L>: right track coordinates
        TruePred        % <2,N>: true state predictions
        AppPred         % <2,N>: approximated state predictions
        ref             % <2,N>: references for each prediction
    end
    
    methods
        function obj = SingleTrackAnimation(track_l, track_r)
            obj.track_l = track_l;
            obj.track_r = track_r;
        end
        
        function initGraphics(obj)
            
        end
        
        function updateGraphics(obj)
        end
    end
end

