classdef singletrackModel
%--------------------------------------------------------------------------
%   xk+1 = fd(xk,uk) + Bd*d(zk),    
%
%       where: zk=Bz*xk and  d~N(mean_d(zk),var_d(zk))
%
%   
%   x = [s, ds, th, dth]'   carriage position and pole angle (and derivatives)
%   u = [F]'                force on the carriage and torque on the pole joint
%   
%--------------------------------------------------------------------------
    
    properties
        cf
        cr
        m   = 1239; % vehicle mass
        Iz  = 1752; % vehicle moment of inertia (yaw axis)
        g   = 9.81; % gravitation
        lf  = 1.19016; % distance of the front wheel to the center of mass 
        lr  = 1.37484; % distance of the rear wheel to the center of mass
        i_g = [3.91 2.002 1.33 1 0.805]; % transmissions of the 1st ... 5th gear

        d
    end
    
    properties(SetAccess=private)
        Bd = [zeros(6,4);
              eye(4)]       % xk+1 = fd(xk,uk) + Bd*d(zk)
        Bz = eye(10)        % z = Bz*x     
        n = 10              % number of outputs x(t)
        m = 5               % number of inputs u(t)
        nd                  % output dimension of du(z)
        nz                  % dimension of z(t)
    end
    
    methods
        function obj = singletrackModel()
            % discretize model
        end
        
        function [mu_xdot, var_xdot] = f(obj, x, u)
            
            sx = x(1);
            sy = x(2);
            v = x(3);
            beta = x(4);
            psi = x(5);
            omega = x(6);
            x_dot = x(7);
            y_dot = x(8);
            psi_dot = x(9);
            varphi_dot = x(10);
            
            delta = u(1);  % steering angle
            G     = u(2);  % gear
            Fb    = u(3);  % brake force
            zeta  = u(4);  % brake force distribution
            phi   = u(5);  % acc pedal position
            
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
            
            
            % calculate xdot mean and covariance
            mu_xdot  = [x_dot; y_dot; v_dot; beta_dot; psi_dot; omega_dot; x_dot_dot; y_dot_dot; psi_dot_dot; varphi_dot_dot];
            var_xdot = zeros(obj.n);
        end
        
        
        function [mu_xkp1, var_xkp1] = fd (obj,xk,uk,dt)
        %------------------------------------------------------------------
        %   Discrete time dynamics of the inverted pendulum (including
        %   disturbance)
        %------------------------------------------------------------------
            % calculate continous time dynamics
            [mu_xdot, var_xdot] = obj.f(xk,uk);
            
            % discretize mean and variance
            mu_xkp1  = xk + dt * mu_xdot;
            var_xkp1 =      dt * var_xdot;
        end
        
        
        function [mu_xkp1,var_xkp1] = xkp1 (obj,xk,uk,dt)
        %------------------------------------------------------------------
        %   Discrete time dynamics of the inverted pendulum (including
        %   disturbance)
        %
        %       xk+1 = fd(xk,uk) + Bd*d(zk)
        %
        %------------------------------------------------------------------
            % calculate continous time dynamics
            [mu_fd, var_fd] = obj.fd(xk,uk,dt);
            
            % evaluate disturbance
            z = obj.Bz * xk;
            [mu_d, var_d] = obj.d(z);
            
            % discretize mean and variance
            mu_xkp1  = mu_fd  + obj.Bd * mu_d;
            var_xkp1 = var_fd + obj.Bd * var_d * obj.Bd';
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

