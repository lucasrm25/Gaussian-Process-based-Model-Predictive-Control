%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -
%------------------------------------------------------------------

% % SAVE CODE FOR LATER
% 
% i_g  = [3.91 2.002 1.33 1 0.805] % transmissions of the 1st ... 5th gear
% i_0  = 3.91         % motor transmission
% R    = 0.302        % wheel radius
% nmax = 4800*2*pi/60 % maximum motor rotation
% 
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



classdef MotionModelGP_SingleTrack < MotionModelGP
%--------------------------------------------------------------------------
%   xk+1 = fd(xk,uk) + Bd * ( d(zk) + w ),    
%
%       where: zk = Bz*xk,
%              d ~ N(mean_d(zk),var_d(zk))
%              w ~ N(0,sigmaw)
%
%   
%   x = [x          (x position), 
%        y          (y position), 
%        v          (velocity), 
%        beta       (side slip angle), 
%        psi        (yaw angle), 
%        omega      (yaw rate), 
%        x_dot      (longitudinal velocity), 
%        y_dot      (lateral velocity), 
%        psi_dot    (yaw rate (redundant)), 
%        varphi_dot (wheel rotary frequency)]'   
%
%   u = [delta      (steering angle), 
%        G          (gear), 
%        F_b        (brake force), 
%        zeta       (brake force distribution), 
%        phi        (acc pedal position)]'               
%   
%--------------------------------------------------------------------------
 
    properties
        M=200; % vehicle mass
        g=9.81; % gravitation
        l_f=1.19016; % distance of the front wheel to the center of mass 
        l_r=1.37484; % distance of the rear wheel to the center of mass
        %l=l_f+l_r; % vehicle length (obsolete)
        R=0.302; % wheel radius
        I_z=500; % vehicle moment of inertia (yaw axis)
        I_R=1.5; % wheel moment of inertia
        i_g=[3.91 2.002 1.33 1 0.805]; % transmissions of the 1st ... 5th gear
        i_0=3.91; % motor transmission
        B_f=10.96; % stiffnes factor (Pacejka) (front wheel)
        C_f=1.3; % shape factor (Pacejka) (front wheel)
        D_f=4560.4; % peak value (Pacejka) (front wheel)
        E_f=-0.5; % curvature factor (Pacejka) (front wheel)
        B_r=12.67; %stiffnes factor (Pacejka) (rear wheel)
        C_r=1.3; %shape factor (Pacejka) (rear wheel)
        D_r=3947.81; %peak value (Pacejka) (rear wheel)
        E_r=-0.5; % curvature factor (Pacejka) (rear wheel)
        f_r_0=0.009; % coefficient (friction)
        f_r_1=0.002; % coefficient (friction)
        f_r_4=0.0003; % coefficient (friction)
    end
    
    properties(SetAccess=private)
        Bd = [0 0 0 0 0 0 0 0 0 0]';% xk+1 = fd(xk,uk) + Bd*d(zk)
        Bz = zeros(10)          % z = Bz*x     
        n = 10                      % number of outputs x(t)
        m = 5                       % number of inputs u(t)
    end
    
    methods
        function obj = MotionModelGP_SingleTrack(d,sigmaw)
        %------------------------------------------------------------------
        %   object constructor
        %------------------------------------------------------------------
            % call superclass constructor
            obj = obj@MotionModelGP(d,sigmaw);
        end
        
        function x = clip(~,x,lb,ub)
            x = min(max(x,lb),ub);
        end
        
        function u = clipInputs(obj,u)
            u(1) = obj.clip( u(1), -deg2rad(30), deg2rad(30) );      % steering angle
            u(2) = round(obj.clip( u(2), 1, 5) );   % gear
            u(3) = obj.clip( u(3), 0, 5000);       % brake force
            u(4) = obj.clip( u(4), 0, 1);           % brake force distribution
            u(5) = obj.clip( u(5), 0, 1);           % acc pedal position
        end
        
        function [xdot, grad_xdot] = f (obj, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics of the single track (including
        %   disturbance):
        %------------------------------------------------------------------
            
            %--------------------------------------------------------------
            % Inputs
            %--------------------------------------------------------------
            u = clipInputs(obj,u);
            
            delta = u(1); % steering angle 
            G     = u(2); % gear 1 ... 5
            F_b   = u(3); %braking force
            zeta  = u(4); % braking force distribution
            phi   = u(5); % gas pedal position

            %--------------------------------------------------------------
            % State Vector
            %--------------------------------------------------------------
            %x = x(1); % x position (obsolete)
            %y = x(2); % y position (obsolete)
            v = x(3); % velocity
            beta = x(4); % side slip angle
            psi = x(5); % yaw angle
            omega = x(6); % yaw rate
            %x_dot=x(7); % longitudinal velocity (obsolete)
            %y_dot=x(8); % lateral velocity (obsolete)
            psi_dot = x(9); % yaw rate (redundant)
            varphi_dot = x(10); % wheel rotary frequency
            
            %--------------------------------------------------------------
            % Dynamics
            %--------------------------------------------------------------

            %--------------------------------------------------------------
            % Slip
            %--------------------------------------------------------------
            
            if v < 3             % do not allow sideslip when velocity is low
                beta = 0;
            end

            %slip angles and steering
            a_f=delta-atan((obj.l_f*psi_dot-v*sin(beta))/(v*cos(beta))); % front slip angle
            a_r=atan((obj.l_r*psi_dot+v*sin(beta))/(v*cos(beta))); %rear slip angle

            
            if isnan(a_f) % front slip angle well-defined?
                a_f=0; % recover front slip angle
            end
            if isnan(a_r) % rear slip angle well-defined
                a_r=0; % recover rear slip angle
            end
            %wheel slip
            % wheel can only rotate forward
            if varphi_dot < 0
                varphi_dot = 0;
            end
            if v<=obj.R*varphi_dot % traction slip? (else: braking slip)
                S=1-(v/(obj.R*varphi_dot)); %traction wheel slip
            else
                S=1-((obj.R*varphi_dot)/v); % braking slip
            end
            if isnan(S) % wheel slip well-defined?
                S=0; % recover wheel slip
            end
            % S=0; % neglect wheel slip

            
            %--------------------------------------------------------------
            % traction, friction, braking
            %--------------------------------------------------------------
            n=v*obj.i_g(G)*obj.i_0*(1/(1-S))/obj.R; % motor rotary frequency
            if isnan(n) % rotary frequency well defined?
                n=0; %recover rotary frequency
            end
            if n>(4800*pi)/30 % maximal rotary frequency exceeded?
                n=(4800*pi)/30; % recover maximal rotary frequency
            end
            
            T_M=200*phi*(15-14*phi)-200*phi*(15-14*phi)*(((n*(30/pi))^(5*phi))/(4800^(5*phi))); % motor torque
            M_wheel=obj.i_g(G)*obj.i_0*T_M; % wheel torque
            F_w_r=(obj.M*obj.l_f*obj.g)/(obj.l_f+obj.l_r); % weight rear
            F_w_f=(obj.M*obj.l_r*obj.g)/(obj.l_f+obj.l_r); % weight front
            f_r=obj.f_r_0+obj.f_r_1*(abs(v)*3.6)/100+obj.f_r_4*((abs(v)*3.6)/100)^4; % approximate friction
            F_b_r=zeta*F_b; % braking force rear
            F_b_f=F_b*(1-zeta); % braking force front
            F_f_r=f_r*F_w_r; % friction rear
            F_f_f=f_r*F_w_f; % friction front
            
            F_x_r=(M_wheel/obj.R)-sign(v*cos(beta))*F_b_r-sign(v*cos(beta))*F_f_r; % longitudinal force rear wheel
            F_x_f=-sign(v*cos(beta))*F_b_f-sign(v*cos(beta))*F_f_f; % longitudinal force front wheel
            F_y_r=obj.D_r*sin(obj.C_r*atan(obj.B_r*a_r-obj.E_r*(obj.B_r*a_r -atan(obj.B_r*a_r)))); % rear lateral force
            F_y_f=obj.D_f*sin(obj.C_f*atan(obj.B_f*a_f-obj.E_f*(obj.B_f*a_f -atan(obj.B_f*a_f)))); % front lateral force

            %--------------------------------------------------------------
            % Output
            %--------------------------------------------------------------
            % vector field (right-hand side of differential equation)
            x_dot=v*cos(psi-beta); % longitudinal velocity
            y_dot=v*sin(psi-beta); % lateral velocity
            v_dot=( F_x_r*cos(beta) + F_x_f*cos(delta+beta) - F_y_r*sin(beta) - F_y_f*sin(delta+beta) )/obj.M; % acceleration
                           
            if abs(v) < 3
                beta_dot = omega;   % we can not have huge sideslip angle if vehicle is not moving
            else
                beta_dot = omega-(F_x_r*sin(beta)+F_x_f*sin(delta+beta)+F_y_r*cos(beta) +F_y_f*cos(delta+beta))/(obj.M*v); % side slip rate
            end
            
            psi_dot=omega; % yaw rate
            omega_dot=(F_y_f*obj.l_f*cos(delta)-F_y_r*obj.l_r +F_x_f*obj.l_f*sin(delta))/obj.I_z; % yaw angular acceleration
            x_dot_dot=(F_x_r*cos(psi)+F_x_f*cos(delta+psi)-F_y_f*sin(delta+psi) -F_y_r*sin(psi))/obj.M; % longitudinal acceleration
            y_dot_dot=(F_x_r*sin(psi)+F_x_f*sin(delta+psi)+F_y_f*cos(delta+psi) +F_y_r*cos(psi))/obj.M; % lateral acceleration
            psi_dot_dot=(F_y_f*obj.l_f*cos(delta)-F_y_r*obj.l_r +F_x_f*obj.l_f*sin(delta))/obj.I_z; % yaw angular acceleration
            varphi_dot_dot=(F_x_r*obj.R)/obj.I_R; % wheel rotary acceleration

            if v_dot< 0
                1;
            end
            if varphi_dot + varphi_dot_dot*0.1 < 0
                varphi_dot_dot = -varphi_dot/0.1;
            end

            %--------------------------------------------------------------
            % write outputs
            %--------------------------------------------------------------
            xdot=[x_dot;y_dot;v_dot;beta_dot;psi_dot;omega_dot;x_dot_dot;y_dot_dot;psi_dot_dot;varphi_dot_dot]; % left-hand side
            grad_xdot = zeros(obj.n);
            
            if any(isnan(xdot)) || any(isinf(xdot)) || any(imag(xdot)~=0)
                error('Single Track Model evaluated to Inf of NaN... CHECK MODEL!!!')
            end
        end
    end
end

