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
        [xdot, grad_xdot] = f (obj, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics.
        %   out:
        %       xdot: <n,1> time derivative of x given x and u
        %       grad_xdot: <n,n> gradient of xdot w.r.t. x
        %-----------------------------------------------------------------  
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
        
        function [xkp1, grad_xkp1] = fd (obj, xk, uk, dt)
        %------------------------------------------------------------------
        %   Discrete time dynamics (ODE1 Euler approximation)
        %   args:
        %       xkp1: <n,1> state prediction (without disturbance model)
        %       grad_xkp1: <n,n> gradient of xkp1 w.r.t. xk
        %------------------------------------------------------------------
            % calculate continous time dynamics
            [xdot, grad_xdot] = obj.f(xk,uk);
            
            % discretize
            xkp1      = xk + dt * xdot;
            grad_xkp1 = eye(obj.n) + dt * grad_xdot;
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
            [fd, grad_fd] = obj.fd(mu_xk,uk,dt);
            
            % calculate grad_{x,d,w} xk+1
            grad_xkp1 = [grad_fd; obj.Bd'; obj.Bd'];
            
            % evaluate disturbance
            z = obj.Bz * mu_xk;
            [mu_d, var_d] = obj.d(z);
            [mu_w, var_w] = obj.w(z);
            
            % a) Mean Equivalent Approximation:
            % var_x_d_w = blkdiag(var_xk, var_d, var_w);
            
            % predict mean and variance (Extended Kalman Filter)
            mu_xkp1  = fd  + obj.Bd * ( mu_d + mu_w );
            var_xkp1 = zeros(obj.n); %grad_xkp1' * var_x_d_w * grad_xkp1;
        end
        
    end
end
