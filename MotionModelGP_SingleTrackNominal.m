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
%        track_dist       (distance traveled in the track centerline)
%        ]
%
%   u = [delta          (steering angle),
%        T              (wheel torque gain),  -1=max.braking, 1=max acc.
%        track_vel      (velocity in the track centerline)
%       ]
%   
%--------------------------------------------------------------------------
 
    properties
        %******************************************************************
        % PLEASE REMEMBER TO REGENERATE THE GRADIENT FUNCTIONS EVERY TIME A
        % PARAMETER IS CHANGED !!!! SEE FUNCTION generate_grad_functions()
        %******************************************************************
        M    = 500      % vehicle mass
        I_z  = 800      % vehicle moment of inertia (yaw axis)
        g    = 9.81     % gravitation
        l_f  = 0.9      % distance of the front wheel to the center of mass 
        l_r  = 1.5      % distance of the rear wheel to the center of mass
        
        deltamax    = deg2rad(30)   % maximum steering amplitude
        % deltadotmax = deg2rad(20)   % maximum steering velocity amplitude
        
        maxbrakeWForce = 6000 % = 2*g*M;  % allow ~ 2g brake
        maxmotorWForce = 6000 % = 1*g*M;  % allow ~ 1g acc
        
        c_f = 15000 % = 1*g*M/deltamax  % front coornering stiffness (C*delta=Fy~M*a)
        c_r = 20000 % = 2*g*M/deltamax  % rear coornering stiffness
        %******************************************************************
        % PLEASE REMEMBER TO REGENERATE THE GRADIENT FUNCTIONS EVERY TIME A
        % PARAMETER IS CHANGED !!!! SEE FUNCTION generate_grad_functions()
        %******************************************************************
    end
    
    properties(SetAccess=private)
        Bd = zeros(7,1);   % xk+1 = fd(xk,uk) + Bd*(d(Bz*xk)+w)
        Bz = eye(7)        % z = Bz*x     
        n = 7              % number of outputs x(t)
        m = 3              % number of inputs u(t)
    end
    
    
    methods
        function x = sclip(obj,x,lb,ub)
            % Smooth (differentiable) clip (saturation) function
            x = x.*obj.gez(x-lb).*obj.lez(x-ub) + ub*obj.gez(x-ub) + lb*obj.lez(x-lb);
        end
    end
    methods(Static)
        function x = clip(x,lb,ub)
            % standard nonsmooth clip (saturation) function
            x = min(max(x,lb),ub);
        end
        function x = srec(x,lb,ub)
            % Smooth rectangular function
            alpha = 50; % the larger the sharper the rectangular function
            x = 0.5*(tanh((x-lb)*alpha)-tanh((x-ub)*alpha));
        end
        function x = gez(x)
            % Smooth >=0 boolean function
            alpha = 50; % the larger the sharper the clip function
            x = (1+exp(-alpha*x)).^-1;
        end
        function x = lez(x)
            % Smooth <=0 boolean function
            alpha = 50; % the larger the sharper the clip function
            x = 1-(1+exp(-alpha*x)).^-1;
        end
        function x = ssign(x)
            % Smooth sign(x) boolean function
            alpha = 100; % the larger the sharper the clip function
            x = tanh(alpha*x);
        end
    end
    
    methods
        function obj = MotionModelGP_SingleTrackNominal(d,sigmaw)
        %------------------------------------------------------------------
        %   object constructor. Create model and report model stability
        %   analysis
        %------------------------------------------------------------------
            % call superclass constructor
            obj = obj@MotionModelGP(d,sigmaw);
            %------------------------------------------------------------------
            %   report created single track model dynamics
            %------------------------------------------------------------------
            fprintf('Single track model created!!! Summary:\n');
            % vehicle size
            l = obj.l_f + obj.l_r;
            % Eigenlenkgradient = Yaw gradient (d_delta/d_ay)
            EG = obj.M/l*(obj.l_r/obj.c_f - obj.l_f/obj.c_r);
            fprintf(2,'\tYaw gradient EG=(d_delta/d_ay) ~ %.1f [deg/g]\n',rad2deg(EG)*obj.g);
            if rad2deg(EG)*obj.g < 10
                fprintf(2,'\tBE CAREFULL... EG is too low (< 10 [deg/g])... increase l_r,c_r or decrease l_f,c_f\n');
            end
            if EG < 0
                v_cr = sqrt( obj.c_f*obj.c_r*l^2 / (obj.M*(obj.c_f*obj.l_f - obj.c_r*obj.l_r)));
                fprintf(2,'\tVehicle is intrinsically OVERSTEER (c_r*l_r < c_f*l_f)\n');
                fprintf('\tCritical velocity v_cr ~ %.1f [m/s]\n', v_cr);
                fprintf('\tMax ss. Yaw (at delta_max) at v= 5[m/s] ~ %.1f [deg/s]\n',rad2deg(obj.deltamax)* 5/(l+EG* 5^2));
                fprintf('\tMax ss. Yaw (at delta_max) at v=10[m/s] ~ %.1f [deg/s]\n',rad2deg(obj.deltamax)*10/(l+EG*10^2));
                fprintf('\tMax ss. Yaw (at delta_max) at v=15[m/s] ~ %.1f [deg/s]\n',rad2deg(obj.deltamax)*15/(l+EG*15^2));
            else
                v_ch = sqrt(l/EG);
                fprintf(2,'\tVehicle is intrinsically UNDERSTEER (c_r*l_r > c_f*l_f)\n');
                fprintf('\tCharacteristic velocity v_ch ~ %.1f [m/s]\n',v_ch);
                fprintf('\tMax ss. Yaw rate gain (at v_ch)          (psidot/delta)ss,max ~ %.1f [1/s]\n',v_ch/(l+EG*v_ch^2));
                fprintf('\tMax ss. Yaw rate (at v_ch and delta_max) (psidot)ss,max~ %.1f [deg/s]\n',rad2deg(obj.deltamax)*v_ch/(l+EG*v_ch^2));
            end
            fprintf('\tax_{max} ~ %.1f [g]\n',obj.maxmotorWForce/obj.M/obj.g);
            fprintf('\tax_{min} ~ -%.1f [g]\n',obj.maxbrakeWForce/obj.M/obj.g);
        end
        
        function xdot = f (obj, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics of the single track (including
        %   disturbance):
        %------------------------------------------------------------------

            %--------------------------------------------------------------
            % State Vector
            %--------------------------------------------------------------
            I_x         = x(1);
            I_y         = x(2);
            psi         = x(3);
            V_vx        = x(4);
            V_vy        = x(5);
            psi_dot     = x(6);
            track_dist  = x(7);
            beta = atan2(V_vy,V_vx);
            
            %--------------------------------------------------------------
            % Inputs
            %--------------------------------------------------------------
            delta     = u(1);   % steering angle
            T         = u(2);   % wheel torque gain,  -1=max.braking, 1=max acc.
            track_vel = u(3);   % track centerline velocity
            
            %--------------------------------------------------------------
            % Saturate inputs and 
            %--------------------------------------------------------------
            % saturate steering angle
            % (DEPRECATED - NOT DIFFERENTIABLE) delta = obj.clip(delta, -obj.deltamax, obj.deltamax);
            delta = obj.sclip(delta, -obj.deltamax, obj.deltamax);
            
            % saturate pedal input
            % delta_dot = obj.clip( delta_dot, -obj.deltadotmax, obj.deltadotmax);
            % (DEPRECATED - NOT DIFFERENTIABLE) T = obj.clip( T, -1, 1);
            T = obj.sclip( T, -1, 1);
            
            %--------------------------------------------------------------
            % Wheel slip angles (slip ration not being used for now)
            %--------------------------------------------------------------
            a_r = atan2(V_vy-obj.l_r*psi_dot,V_vx);
            a_f = atan2(V_vy+obj.l_f*psi_dot,V_vx) - delta;
                
            %--------------------------------------------------------------
            % Tyre forces
            %--------------------------------------------------------------
            % desired total wheel torque to be applied
            % (DEPRECATED - NOT DIFFERENTIABLE) totalWForce = T * ( (T>0)*obj.maxmotorWForce+(T<0)*obj.maxbrakeWForce*sign(V_vx) );
            totalWForce = T*(  (obj.gez(T)).*obj.maxmotorWForce ...
                              +(obj.lez(T)).*obj.maxbrakeWForce.*obj.ssign(V_vx));        
            % longitudinal wheel torque distribution
            zeta = 0.5;
            
            % Wheel forces in wheel coordinates (z-axis points down, x-axis front)
            % This means positive W_Fy_f turns vehicle to the right
            W_Fx_r = zeta * totalWForce;
            W_Fx_f = (1-zeta) * totalWForce;
            W_Fy_r = obj.c_r * a_r;   
            W_Fy_f = obj.c_f * a_f;
            
            % Wheel forces in vehicle coordinates (z-axis points up, x-axis front)
            V_Fx_r = W_Fx_r;
            V_Fx_f = W_Fx_f;
            V_Fy_r = - W_Fy_r;
            V_Fy_f = - W_Fy_f;
            
            %--------------------------------------------------------------
            % Calculate state space time derivatives
            %--------------------------------------------------------------
            % vector field (right-hand side of differential equation)
            I_x_dot = V_vx*cos(psi) - V_vy*sin(psi); % longitudinal velocity
            I_y_dot = V_vx*sin(psi) + V_vy*cos(psi); % lateral velocity
            V_vx_dot = 1/obj.M * (V_Fx_r + V_Fx_f*cos(delta) - V_Fy_f*sin(delta) + V_vy*psi_dot);
            V_vy_dot = 1/obj.M * (V_Fy_r + V_Fx_f*sin(delta) + V_Fy_f*cos(delta) - V_vy*psi_dot);
            psi_dot_dot = 1/obj.I_z * (V_Fy_f*obj.l_f*cos(delta) + V_Fx_f*obj.l_f*sin(delta) - V_Fy_r*obj.l_r);
            track_dist_dot = track_vel; % Traveled distance in the track centerline
                                
            %--------------------------------------------------------------
            % write outputs
            %--------------------------------------------------------------
            xdot  = [I_x_dot; I_y_dot; psi_dot; V_vx_dot; V_vy_dot; psi_dot_dot; track_dist_dot];
            
            %--------------------------------------------------------------
            % check model validity
            %--------------------------------------------------------------
            if ~isa(x,'sym')    % sym input used for code generation
                if any(isnan(xdot)) || any(isinf(xdot)) || any(imag(xdot)~=0)
                    error('Single Track Model evaluated to Inf of NaN... CHECK MODEL!!!')
                end
            end
        end
        
        function gradx = gradx_f(~, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics.
        %   out:
        %       gradx: <n,n> gradient of xdot w.r.t. x
        %------------------------------------------------------------------
            % gradx = zeros(obj.n);
            gradx = CODEGEN_singletrack_gradx_f(x,u);
        end
        
        function gradu = gradu_f(~, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics.
        %   out:
        %       gradu: <m,n> gradient of xdot w.r.t. u
        %------------------------------------------------------------------
            % gradu = zeros(obj.m,obj.n);
            gradu = CODEGEN_singletrack_gradu_f(x,u);
        end
        
        function generate_grad_functions(obj)
        %------------------------------------------------------------------
        %   Generate external files for the evaluation of the gradient of 
        %   the continuous time dynamics. (Make use of symbolic toolbox)
        %   Please ensure that your dynamics only contains smooth
        %   diferentiable functions.
        %   
        %   To generate files, simply call:
        %       nomModel = MotionModelGP_SingleTrackNominal(@(z)deal(0,0), 0);
        %       nomModel.generate_grad_functions() 
        %------------------------------------------------------------------
            syms I_x I_y vpsi V_vx V_vy psi_dot track_dist real
            x = [I_x I_y vpsi V_vx V_vy psi_dot track_dist]';
            syms delta T track_vel real
            u = [delta T track_vel]';
            
            xdot = obj.f(x,u);
            gradx = jacobian(xdot,x)';
            gradu = jacobian(xdot,u)';
            % gradx = simplify(expand(gradx));  % does not work. eqs are too complex
            % gradu = simplify(expand(gradu));  % does not work. eqs are too complex

            matlabFunction(gradx,'Vars',{x,u},'File','CODEGEN_singletrack_gradx_f','Optimize',true);
            matlabFunction(gradu,'Vars',{x,u},'File','CODEGEN_singletrack_gradu_fssss','Optimize',true);
            disp('FINISHED! functions CODEGEN_singletrack_gradx_f and CODEGEN_singletrack_gradu_fssss generated!!')
        end
    end
end
