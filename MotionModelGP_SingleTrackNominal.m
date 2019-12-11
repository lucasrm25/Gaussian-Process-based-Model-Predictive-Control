%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -
%------------------------------------------------------------------

classdef MotionModelGP_SingleTrackNominal < MotionModelGP
%--------------------------------------------------------------------------
%   xk+1 = fd(xk,uk) + Bd * ( d(zk) + w ),    
%
%       where: zk = Bz*xk,
%              d ~ N(mean_d(zk),var_d(zk))
%              w ~ N(0,sigmaw)
%
%   
%   x = [I_x              (x position in global coordinates), 
%        I_y              (y position in global coordinates), 
%        psi              (yaw angle),
%        V_vx             (longitudinal velocity in vehicle coordinates)             
%        V_vy             (lateral velocity in vehicle coordinates)
%        psi_dot          (yaw rate),
%        delta            (steering angle)
%        track_dist       (distance traveled in the track centerline)
%        ]
%
%   u = [delta_dot      (steering angle velocity),
%        T              (wheel torque gain),  -1=max.braking, 1=max acc.
%        zeta           (wheel torque longitudinal distribution), 
%        track_vel      (velocity in the track centerline)
%       ]
%   
%--------------------------------------------------------------------------
 
    properties
        M   = 500 % vehicle mass
        I_z  = 1000 % vehicle moment of inertia (yaw axis)
        g   = 9.81 % gravitation
        l_f  = 1.19016 % distance of the front wheel to the center of mass 
        l_r  = 1.37484 % distance of the rear wheel to the center of mass
        i_g = [3.91 2.002 1.33 1 0.805] % transmissions of the 1st ... 5th gear
        i_0 = 3.91 % motor transmission
        R = 0.302 % wheel radius
        nmax = 4800*2*pi/60 % maximum motor rotation
        
        deltamax = deg2rad(30)  % maximum steering amplitude
        deltadotmax = deg2rad(10) % maximum steering velocity amplitude
        
        maxbrakeWForce % = -2*g*M;  % allow ~ 2g brake
        maxmotorWForce % =  1*g*M;  % allow ~ 1g acc
        
        cf  % = 1*g*M/deltamax  % front coornering stiffness (C*delta=Fy~M*a)
        cr  % = 2*g*M/deltamax  % rear coornering stiffness
    end
    
    properties(SetAccess=private)
        Bd = zeros(8,1);   % xk+1 = fd(xk,uk) + Bd*(d(Bz*xk)+w)
        Bz = eye(8)        % z = Bz*x     
        n = 8              % number of outputs x(t)
        m = 4              % number of inputs u(t)
    end
    
    methods
        function obj = MotionModelGP_SingleTrackNominal(d,sigmaw)
        %------------------------------------------------------------------
        %   object constructor
        %------------------------------------------------------------------
            % call superclass constructor
            obj = obj@MotionModelGP(d,sigmaw);
            
            obj.maxbrakeWForce = -2*obj.g*obj.M;  % allow ~ 2g brake
            obj.maxmotorWForce =  1*obj.g*obj.M;  % allow ~ 1g acc

            obj.cf  = 1*obj.g*obj.M/obj.deltamax;  % front coornering stiffness (C*delta=Fy~M*a)
            obj.cr  = 2*obj.g*obj.M/obj.deltamax;  % rear coornering stiffness
        end
        
        function x = clip(~,x,lb,ub)
            x = min(max(x,lb),ub);
        end
        
        function [xdot, grad_xdot] = f (obj, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics of the single track (including
        %   disturbance):
        %------------------------------------------------------------------


            %--------------------------------------------------------------
            % State Vector
            %--------------------------------------------------------------
            I_x     = x(1);
            I_y     = x(2);
            psi     = x(3);
            V_vx    = x(4);
            V_vy    = x(5);
            psi_dot = x(6);
            delta   = x(7);
            track_dist = x(8);
            
            % saturate steering angle
            delta = obj.clip(delta, -obj.deltamax, obj.deltamax);
            
            % calculate sideslip angle
            beta = atan2(V_vy,V_vx);
            
            %--------------------------------------------------------------
            % Inputs
            %--------------------------------------------------------------
            delta_dot = u(1);   % steering angle velocity
            T         = u(2);   % wheel torque gain,  -1=max.braking, 1=max acc.
            zeta      = u(3);   % wheel torque longitudinal distribution
            track_vel = u(4);   % track centerline velocity
            
            % saturate inputs to valid ranges
            delta_dot = obj.clip( delta_dot, -obj.deltadotmax, obj.deltadotmax);
            T         = obj.clip( T, -1, 1);
            zeta      = obj.clip( zeta,  0, 1);
            
            %--------------------------------------------------------------
            % Traveled distance in the track centerline
            %--------------------------------------------------------------
            track_dist_dot = track_vel;
            
            %--------------------------------------------------------------
            % Slip
            %--------------------------------------------------------------
            % slip angles and steering
            if V_vx > 3
                a_r = atan2(V_vy-obj.l_r*psi_dot,V_vx);
                a_f = atan2(V_vy+obj.l_f*psi_dot,V_vx) - delta;
            else
                a_r = 0;
                a_f = - delta;
            end
                
            %--------------------------------------------------------------
            % Tyre forces
            %--------------------------------------------------------------
            % % motor rotary frequency
            % n = V_vx/obj.R * obj.i_g(G) * obj.i_0; 
            % if n > 0 && n < obj.nmax
                % % motor torque
                % T_M = 200*phi*(15-14*phi)-200*phi*(15-14*phi)*(((n*(30/pi))^(5*phi))/(4800^(5*phi)));
            %else
                %T_M = 0;    % motor outside rotation range
            %end
            % % wheel torque
            % T_W = T_M * obj.i_g(G) * obj.i_0;
            % W_Fx_r =     -zeta*F_b*(sign(V_vx)) + T_W/obj.R;  
            % W_Fx_f = -(1-zeta)*F_b*(sign(V_vx));
            % W_Fy_r = obj.cr * a_r;
            % W_Fy_f = obj.cf * a_f;
            % W_Fx_r = obj.clip(W_Fx_r,-5000,5000);
            % W_Fx_f = obj.clip(W_Fx_f,-5000,5000);
            % W_Fy_r = obj.clip(W_Fy_r,-5000,5000);
            % W_Fy_f = obj.clip(W_Fy_f,-5000,5000);
            
            % desired total wheel torque to be applied
            totalWForce = T * ( (T>0)*obj.maxmotorWForce+(T<0)*obj.maxbrakeWForce*sign(V_vx) );
            
            % wheel forces in wheel coordinates
            W_Fx_r = zeta * totalWForce;
            W_Fx_f = (1-zeta) * totalWForce;
            W_Fy_r = obj.cr * a_r;
            W_Fy_f = obj.cf * a_f;
            
            %--------------------------------------------------------------
            % Output
            %--------------------------------------------------------------
            % vector field (right-hand side of differential equation)
            I_x_dot = V_vx*cos(psi) - V_vy*sin(psi); % longitudinal velocity
            I_y_dot = V_vx*sin(psi) + V_vy*cos(psi); % lateral velocity
            
            V_vx_dot = 1/obj.M * (W_Fx_r + W_Fx_f*cos(delta) - W_Fy_f*sin(delta)) + V_vy*psi_dot;
            V_vy_dot = 1/obj.M * (W_Fy_r + W_Fx_f*sin(delta) + W_Fy_f*cos(delta)) - V_vy*psi_dot;
            
             % yaw angular acceleration;
            psi_dot_dot = 1/obj.I_z * (W_Fy_f*obj.l_f*cos(delta) + W_Fx_f*obj.l_f*sin(delta) - W_Fy_r*obj.l_r);
                    
            %--------------------------------------------------------------
            % write outputs
            %--------------------------------------------------------------
            xdot  = [I_x_dot; I_y_dot; psi_dot; V_vx_dot; V_vy_dot; psi_dot_dot; delta_dot; track_dist_dot];
            grad_xdot = zeros(obj.n);
            
            if any(isnan(xdot)) || any(isinf(xdot)) || any(imag(xdot)~=0)
                error('Single Track Model evaluated to Inf of NaN... CHECK MODEL!!!')
            end
        end
    end
end

