classdef singletracknominal
    %SINGLETRACKNOMINAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cf
        cr
        m   = 1239; % vehicle mass
        Iz  = 1752; % vehicle moment of inertia (yaw axis)
        g   = 9.81; % gravitation
        lf  = 1.19016; % distance of the front wheel to the center of mass 
        lr  = 1.37484; % distance of the rear wheel to the center of mass
        i_g = [3.91 2.002 1.33 1 0.805]; % transmissions of the 1st ... 5th gear
        deltamax
        vmax
    end
    
    methods
        function obj = singletracknominal()
            % discretize model
        end
        
        function xkp1 = f(obj,xk)
            xkp1 = xk;
        end
    end
end

