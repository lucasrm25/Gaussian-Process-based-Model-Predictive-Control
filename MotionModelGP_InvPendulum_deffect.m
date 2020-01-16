%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -
%------------------------------------------------------------------


classdef MotionModelGP_InvPendulum_deffect < MotionModelGP_InvPendulum_nominal
%--------------------------------------------------------------------------
%   Inverted pendulum dynamics based on the nominal model but with some 
%   deffect in the actuators
%
%   xk+1 = fd_deffect(xk,uk) + Bd * ( d(zk) + w ),    
%
%       where: zk = [Bz_x*xk ; Bz_u*uk],
%              d ~ N(mean_d(zk),var_d(zk))
%              w ~ N(0,sigmaw)
%
%   
%   x = [s, ds, th, dth]'   carriage position and pole angle (and derivatives)
%   u = [F]'                force on the carriage and torque on the pole joint
%   
%--------------------------------------------------------------------------

    methods
        function obj = MotionModelGP_InvPendulum_deffect(Mc, Mp, b, I, l, d, sigmaw)
            obj = obj @ MotionModelGP_InvPendulum_nominal(Mc, Mp, b, I, l, d, sigmaw);
        end
        
        
        function xdot = f (obj, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics of the inverted pendulum but with some 
        %   deffect (additional disturbance)
        %   args:
        %       x: <n,1>
        %       u: <m,1>
        %   out:
        %       xdot: <n,1> time derivative of x given x and u
        %------------------------------------------------------------------
            % get dynamics from nominal model
            xdot = f @ MotionModelGP_InvPendulum_nominal(obj,x,u);
            
            % add deffect
            xdot(3) = xdot(3) + (0.1 * x(3) - 0.01*x(4) + deg2rad(3)) *10;
            
            % xdot(2) = xdot(2) + ( -0.1 * x(1) + deg2rad(3));
        end
        
%         function [xkp1, gradx_xkp1] = fd (obj, xk, uk, dt)
%             % get dynamics from nominal model
%             [xkp1, gradx_xkp1] = fd @ MotionModelGP_InvPendulum_nominal(obj, xk, uk, dt);
%             % deffect
%             xkp1(3) = xkp1(3) + (0.1 * xk(3) - 0.01*xk(4) + deg2rad(3)) *1;
%         end
        
    end

end