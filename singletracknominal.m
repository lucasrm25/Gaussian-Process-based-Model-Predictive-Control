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
        
        function xkp1 = f(obj, xk, uk, dt)
            
            x = xk(1);
            y = xk(2);
            v = xk(3);
            beta = xk(4);
            psi = xk(5);
            omega = xk(6);
            x_dot = xk(7);
            y_dot = xk(8);
            psi_dot = xk(9);
            varphi_dot = xk(10);
            
            delta = uk(1);  % steering angle
            G     = uk(2);  % gear
            Fb    = uk(3);  % brake force
            zeta  = uk(4);  % brake force distribution
            phi   = uk(5);  % acc pedal position
            
            T_M=200*phi*(15-14*phi)-200*phi*(15-14*phi)*(((n*(30/pi))^(5*phi))/(4800^(5*phi))); % motor torque
            
            F_x_r = T_M - F_b;  % longitudinal force rear wheel
            F_x_f = 0;          % longitudinal force front wheel
            F_y_r = Cr * a_r;   % rear lateral force
            F_y_f = Cf * a_f;   % front lateral force

            % vector field (right-hand side of differential equation)
            
            x_dot = v*cos(psi-beta); % longitudinal velocity
            y_dot = v*sin(psi-beta); % lateral velocity
            v_dot = (F_x_r*cos(beta)+F_x_f*cos(delta+beta)-F_y_r*sin(beta) -F_y_f*sin(delta+beta))/obj.m; % acceleration
            beta_dot = omega-(F_x_r*sin(beta)+F_x_f*sin(delta+beta)+F_y_r*cos(beta) +F_y_f*cos(delta+beta))/(obj.m*v); % side slip rate
            psi_dot  = omega; % yaw rate
            omega_dot=(F_y_f*l_f*cos(delta)-F_y_r*l_r +F_x_f*l_f*sin(delta))/obj.Iz; % yaw angular acceleration
            x_dot_dot=(F_x_r*cos(psi)+F_x_f*cos(delta+psi)-F_y_f*sin(delta+psi) -F_y_r*sin(psi))/obj.m; % longitudinal acceleration
            y_dot_dot=(F_x_r*sin(psi)+F_x_f*sin(delta+psi)+F_y_f*cos(delta+psi) +F_y_r*cos(psi))/obj.m; % lateral acceleration
            psi_dot_dot=(F_y_f*l_f*cos(delta)-F_y_r*l_r  +F_x_f*l_f*sin(delta))/obj.Iz; % yaw angular acceleration
            varphi_dot_dot=(F_x_r*R)/obj.Iz; % wheel rotary acceleration
            
            xdot = [x_dot; y_dot; v_dot; beta_dot; psi_dot; omega_dot; x_dot_dot; y_dot_dot; psi_dot_dot; varphi_dot_dot];
            
            % discretize continuous time model
            xkp1 = xk + dt * xdot;
        end
        
        
        function r = ref(obj, tk, xk, t_r, t_l)
            %     xk = [2,1]';
            % calculate trajectory center line
            t_c = (t_r + t_l)/2;
            % find closest trajectory point w.r.t. the vehicle
            [~,idx] = min( pdist2(xk',t_c,'seuclidean',[1 1].^0.5).^2 );
            % set target as 3 poins ahead
            idx_target = idx +3;
            % loop around when track is over
            idx_target = mod(idx_target, size(t_c,1));
            % return reference signal
            r = t_c(idx_target,:);
        end


    end
end

