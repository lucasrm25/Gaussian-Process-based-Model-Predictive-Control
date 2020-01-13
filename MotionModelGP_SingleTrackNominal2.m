%--------------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -
%--------------------------------------------------------------------------

classdef MotionModelGP_SingleTrackNominal2 < MotionModelGP
%--------------------------------------------------------------------------
%   xk+1 = fd(xk,uk) + Bd * ( d(zk) + w ),    
%
%       where: zk = [Bz_x*xk ; Bz_u*uk],
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
 
    properties(Constant)
        M    = 500      % vehicle mass
        I_z  = 600      % vehicle moment of inertia (yaw axis)
        g    = 9.81     % gravitation
        l_f  = 0.9      % distance of the front wheel to the center of mass 
        l_r  = 1.5      % distance of the rear wheel to the center of mass
        
        deltamax    = deg2rad(30)   % maximum steering amplitude
        % deltadotmax = deg2rad(20)   % maximum steering velocity amplitude
        
        maxbrakeWForce = 8000 % = 2*g*M;  % allow ~ 2g brake
        maxmotorWForce = 4000 % = 1*g*M;  % allow ~ 1g acc
        
        
        % Pacejka lateral dynamics parameters
        B_f = 10.96;              % stiffnes factor (Pacejka) (front wheel)
        C_f = 1.3;                % shape factor (Pacejka) (front wheel)
        D_f = 4560.4;             % peak value (Pacejka) (front wheel)
        E_f = -0.5;               % curvature factor (Pacejka) (front wheel)

        B_r = 12.67;              % stiffnes factor (Pacejka) (rear wheel)
        C_r = 1.3;                % shape factor (Pacejka) (rear wheel)
        D_r = 3947.81;            % peak value (Pacejka) (rear wheel)
        E_r = -0.5;               % curvature factor (Pacejka) (rear wheel)

        
        % keep in mind the dimensions:  xk+1 = fd(xk,uk) + Bd*(d(z)+w)),
        % where z = [Bz_x*x;Bz_u*u] 
        Bz_x = [0 0 0  1 0 0  0;
                0 0 0  0 1 0  0;
                0 0 0  0 0 1  0] 
        Bz_u = [1 0 0;
                0 1 0] 
        Bd = [zeros(3);
              eye(3); 
              0 0 0];               
        n  = 7   % number of outputs x(t)
        m  = 3   % number of inputs u(t)
        nz = 5   % dimension of z(t)
        nd = 3   % output dimension of d(z)
    end
    
    methods
        function obj = MotionModelGP_SingleTrackNominal2(d,sigmaw)
        %------------------------------------------------------------------
        %   object constructor. Create model and report model stability
        %   analysis
        %------------------------------------------------------------------
            % call superclass constructor
            obj = obj@MotionModelGP(d,sigmaw);
            % report single track dynamics analysis
            % obj.analyseSingleTrack();
            
            % add folder CODEGEN to path. Here there will be some functions
            % generated with the method generate_grad_functions()
            addpath(fullfile(pwd,'CODEGEN'))
        end
        
        function xdot = f (obj, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics of the single track (including
        %   disturbance):
        %------------------------------------------------------------------

            %--------------------------------------------------------------
            % Model parameters
            %--------------------------------------------------------------
            g = obj.g;
            % switch to symbolic variables if input is symbolic
            % (to be used with automatic gradient code generation)
            if isa(x,'sym')
                syms M I_z l_f l_r deltamax maxbrakeWForce maxmotorWForce real
                syms B_f C_f D_f E_f real
                syms B_r C_r D_r E_r real
            else
                M = obj.M;
                I_z  = obj.I_z;
                l_f  = obj.l_f;
                l_r  = obj.l_r;
                deltamax = obj.deltamax;
                maxbrakeWForce = obj.maxbrakeWForce;
                maxmotorWForce = obj.maxmotorWForce;
                B_f=obj.B_f; C_f=obj.C_f; D_f=obj.D_f; E_f=obj.E_f;
                B_r=obj.B_r; C_r=obj.C_r; D_r=obj.D_r; E_r=obj.E_r;
            end
        
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
            if isa(x,'sym') % switch to differentiable function
                delta = obj.sclip(delta, -deltamax, deltamax);
            else
                delta = obj.clip(delta, -deltamax, deltamax); % (NOT DIFFERENTIABLE)
            end
            
            % saturate pedal input
            if isa(x,'sym') % switch to differentiable function
                T = obj.sclip( T, -1, 1);
            else
                T = obj.clip( T, -1, 1);    % % (NOT DIFFERENTIABLE) 
            end
            
            %--------------------------------------------------------------
            % Wheel slip angles (slip ration not being used for now)
            %--------------------------------------------------------------
            a_r = atan2(V_vy-l_r*psi_dot,V_vx);
            a_f = atan2(V_vy+l_f*psi_dot,V_vx) - delta;
                
            %--------------------------------------------------------------
            % Tyre forces
            %--------------------------------------------------------------
            % desired total wheel torque to be applied
            if isa(x,'sym') % switch to differentiable function
                totalWForce = T*(  (obj.gez(T)).*maxmotorWForce ...
                              +(obj.lez(T)).*maxbrakeWForce.*obj.ssign(V_vx)); 
            else
                totalWForce = T * ( (T>0)*maxmotorWForce+(T<0)*maxbrakeWForce*sign(V_vx) ); % % (NOT DIFFERENTIABLE) 
            end   
            % longitudinal wheel torque distribution
            zeta = 0.5;
            
            % Wheel forces in wheel coordinates (z-axis points down, x-axis front)
            % This means positive W_Fy_f turns vehicle to the right
            W_Fx_r = zeta * totalWForce;
            W_Fx_f = (1-zeta) * totalWForce;
            
            W_Fy_r = D_r*sin(C_r*atan(B_r*a_r-E_r*(B_r*a_r -atan(B_r*a_r)))); % rear lateral force
            W_Fy_f = D_f*sin(C_f*atan(B_f*a_f-E_f*(B_f*a_f -atan(B_f*a_f)))); % front lateral force
            
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
            V_vx_dot = 1/M * (V_Fx_r + V_Fx_f*cos(delta) - V_Fy_f*sin(delta) + V_vy*psi_dot);
            V_vy_dot = 1/M * (V_Fy_r + V_Fx_f*sin(delta) + V_Fy_f*cos(delta) - V_vy*psi_dot);
            psi_dot_dot = 1/I_z * (V_Fy_f*l_f*cos(delta) + V_Fx_f*l_f*sin(delta) - V_Fy_r*l_r);
            track_dist_dot = track_vel; % Traveled distance in the track centerline
                                
            %--------------------------------------------------------------
            % write outputs
            %--------------------------------------------------------------
            xdot  = [I_x_dot; I_y_dot; psi_dot; V_vx_dot; V_vy_dot; psi_dot_dot; track_dist_dot];
            
            %--------------------------------------------------------------
            % check model integrity
            %--------------------------------------------------------------
            if ~isa(x,'sym')    % sym input used for code generation
                if any(isnan(xdot)) || any(isinf(xdot)) || any(imag(xdot)~=0)
                    error('Single Track Model evaluated to Inf of NaN... CHECK MODEL!!!')
                end
            end
        end
        
        
        function gradx = gradx_f(obj, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics.
        %   out:
        %       gradx: <n,n> gradient of xdot w.r.t. x
        %------------------------------------------------------------------
            % gradx = zeros(obj.n);
            params = [obj.M obj.I_z obj.l_f obj.l_r obj.deltamax ...
                      obj.maxbrakeWForce obj.maxmotorWForce ...
                      obj.B_f obj.C_f obj.D_f obj.E_f ...
                      obj.B_r obj.C_r obj.D_r obj.E_r ]';
            gradx = singletrack_gradx_f(x,u,params);
        end
        
        
        function gradu = gradu_f(obj, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics.
        %   out:
        %       gradu: <m,n> gradient of xdot w.r.t. u
        %------------------------------------------------------------------
            % gradu = zeros(obj.m,obj.n);
            params = [obj.M obj.I_z obj.l_f obj.l_r obj.deltamax ...
                      obj.maxbrakeWForce obj.maxmotorWForce ...
                      obj.B_f obj.C_f obj.D_f obj.E_f ...
                      obj.B_r obj.C_r obj.D_r obj.E_r ]';
            gradu = singletrack_gradu_f(x,u,params);
        end
        
        
        function generate_grad_functions(obj)
        %------------------------------------------------------------------
        %   Generate external files for the evaluation of the gradient of 
        %   the continuous time dynamics. (Make use of symbolic toolbox)
        %   Please ensure that your dynamics only contains smooth
        %   diferentiable functions.
        %   
        %   To generate files, simply call:
        %       nomModel = MotionModelGP_SingleTrackNominal2(@(z)deal(zeros(3,1),zeros(3)), zeros(3));
        %       nomModel.generate_grad_functions() 
        %------------------------------------------------------------------
            syms I_x I_y vpsi V_vx V_vy psi_dot track_dist real
            x = [I_x I_y vpsi V_vx V_vy psi_dot track_dist]';
            syms delta T track_vel real
            u = [delta T track_vel]';
            syms M I_z l_f l_r deltamax maxbrakeWForce maxmotorWForce real
            syms B_f C_f D_f E_f real
            syms B_r C_r D_r E_r real
            
            params = [M I_z l_f l_r deltamax maxbrakeWForce maxmotorWForce ...
                      B_f C_f D_f E_f ...
                      B_r C_r D_r E_r]';
            
            xdot = obj.f(x,u);
            gradx = jacobian(xdot,x)';
            gradu = jacobian(xdot,u)';
            % gradx = simplify(expand(gradx));  % does not work. eqs are too complex
            % gradu = simplify(expand(gradu));  % does not work. eqs are too complex

            % Save model
            folder = fullfile(pwd,'CODEGEN');
            if ~exist(folder,'dir')
                mkdir(folder); 
            end
            addpath(folder)
            
            matlabFunction(gradx,'Vars',{x,u,params},'File',fullfile('CODEGEN','singletrack_gradx_f'),'Optimize',true);
            matlabFunction(gradu,'Vars',{x,u,params},'File',fullfile('CODEGEN','singletrack_gradu_f'),'Optimize',true);
            
            disp('FINISHED! functions CODEGEN_singletrack_gradx_f and CODEGEN_singletrack_gradu_f generated!!')
        end
        
        
        function analyseSingleTrack(obj)
        %------------------------------------------------------------------
        %   report single track model dynamics analysis
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
    end
    
    
end
