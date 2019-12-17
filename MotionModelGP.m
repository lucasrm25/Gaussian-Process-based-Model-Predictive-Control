%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -
%------------------------------------------------------------------

classdef (Abstract) MotionModelGP < handle
%--------------------------------------------------------------------------
%   Abstract class for implementing Motion Model of Plant
%   Intended to be used with GP and NMPC classes
%
%   Please inherit this class and implement all the (Abstract) methods and
%   variables
%
%   xk+1 = fd(xk,uk) + Bd * ( d(zk) + w ),    
%
%       where: zk = Bz*xk,
%              d ~ N(mean_d(zk),var_d(zk))
%              w ~ N(0,sigmaw)
%   
%--------------------------------------------------------------------------

    properties (Abstract, SetAccess=private)
        Bd  % <n,xx> xk+1 = fd(xk,uk) + Bd*d(zk)
        Bz  % <yy,n> z = Bz*x   
        n   % <1>    number of outputs x(t)
        m   % <1>    number of inputs u(t)
    end
    
    properties (SetAccess=private)
        d  % [E[d(z)] , Var[d(z)]] = d(z): disturbace model
        w  % [E[w] , Var[w]] = w: measurement noise
        % this properties are obtained from Set methods
        nd  % output dimension of du(z)
        nz  % dimension of z(t)
    end
    
    methods (Abstract)
        xdot = f (obj, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics.
        %   out:
        %       xdot: <n,1> time derivative of x given x and u
        %------------------------------------------------------------------
        
        gradx = gradx_f(obj, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics.
        %   out:
        %       gradx: <n,n> gradient of xdot w.r.t. x
        %------------------------------------------------------------------
        
        gradu = gradu_f(obj, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics.
        %   out:
        %       gradu: <m,n> gradient of xdot w.r.t. u
        %------------------------------------------------------------------
    end
    
    methods
        function obj = MotionModelGP (d, sigmaw)
        %------------------------------------------------------------------
        %   object constructor
        %   args:
        %       d: evaluates nonlinear motion model mean and covariance 
        %          function [mean_d, var_d] = d(z),   with z = Bz*x
        %       sigmaw: <1> measurement noise covariance
        %------------------------------------------------------------------
            obj.d = d;
            obj.w = @(z) deal(zeros(obj.nd,1),eye(obj.nd)*sigmaw);

            assert(size(obj.Bd,1) == obj.n, 'obj.Bd matrix has wrong number of rows')
            assert(size(obj.Bz,2) == obj.n, 'obj.Bz matrix has wrong number of rows')
        end        
        
        function nd = get.nd(obj)
        %------------------------------------------------------------------
        %   output dimension of d(z)
        %------------------------------------------------------------------
            nd = size(obj.Bd,2);
        end
        
        function nz = get.nz(obj)
        %------------------------------------------------------------------
        %   dimension of zk=Bz*xk
        %------------------------------------------------------------------
            nz = size(obj.Bz,1);
        end
        
        function [xkp1, gradx_xkp1] = fd (obj, xk, uk, dt)
        %------------------------------------------------------------------
        %   Discrete time dynamics (ODE1 Euler approximation)
        %   args:
        %       xkp1: <n,1> state prediction (without disturbance model)
        %       grad_xkp1: <n,n> gradient of xkp1 w.r.t. xk
        %------------------------------------------------------------------
            
            solver = 'ode1';
            
            if strcmp(solver,'ode1')
                %-----------------
                % Fixed step ODE1
                %-----------------
                % calculate continous time dynamics
                xkp1 = xk + dt * obj.f(xk,uk);
                
            elseif strcmp(solver,'ode2')
                %-----------------
                % Fixed step ODE2 (developed by myself)
                %-----------------
                [~,xkp1] = ODE.ode2(@(t,x) obj.f(x,uk), xk, dt, dt);
                
            elseif strcmp(solver,'ode4')
                %-----------------
                % Fixed step ODE4 (developed by myself)
                %-----------------
                [~,xkp1] = ODE.ode4(@(t,x) obj.f(x,uk), xk, dt, dt);
                
            elseif strcmp(solver,'ode23')
                %-----------------
                % Variable step ODE23
                %-----------------
                opts_1 = odeset('Stats','off','RelTol',1e-1,'AbsTol',1e-1);
                [~,xkp1_23] = ode23( @(t,x) obj.f(x,uk), [0 dt], xk, opts_1);
                xkp1 = xkp1_23(end,:)';
            
            elseif strcmp(solver,'ode23')
                %-----------------
                % Variable step ODE23
                %-----------------
                opts_1 = odeset('Stats','off','RelTol',1e-1,'AbsTol',1e-1);
                [~,xkp1_23] = ode23( @(t,x) obj.f(x,uk), [0 dt], xk, opts_1);
                xkp1 = xkp1_23(end,:)';
                
            else
                error('Not implemented');
            end
            
            % for now, gradient is being discretized using a simple ode1
            gradx_xdot = obj.gradx_f(xk,uk);
            gradx_xkp1 = eye(obj.n) + dt * gradx_xdot;
        end
        
        function [mu_xkp1,var_xkp1] = xkp1 (obj, mu_xk, var_xk, uk, dt)
        %------------------------------------------------------------------
        %   State prediction (motion model) using Extended Kalman Filter 
        %   approach
        %
        %       xk+1 = fd(xk,uk) + Bd * ( d(zk) + w ),  zk=Bz*xk
        %
        %------------------------------------------------------------------
            % calculate discrete time dynamics
            [fd, gradx_fd] = obj.fd(mu_xk,uk,dt);
            
            % calculate grad_{x,d,w} xk+1
            grad_xkp1 = [gradx_fd; obj.Bd'; obj.Bd'];
            
            % evaluate disturbance
            z = obj.Bz * mu_xk;
            [mu_d, var_d] = obj.d(z);
            [mu_w, var_w] = obj.w(z);
            
            % a) Mean Equivalent Approximation:
            var_x_d_w = blkdiag(var_xk, var_d, var_w);
            
            % predict mean and variance (Extended Kalman Filter)
            mu_xkp1  = fd  + obj.Bd * ( mu_d + mu_w );
            var_xkp1 = grad_xkp1' * var_x_d_w * grad_xkp1; % zeros(obj.n);
        end
        
    end
end
