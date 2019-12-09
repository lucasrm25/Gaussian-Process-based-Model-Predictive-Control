classdef InvertedPendulumModel < handle
%--------------------------------------------------------------------------
%   xk+1 = fd(xk,uk) + Bd * ( d(zk) + w ),    
%
%       where: zk = Bz*xk,
%              d ~ N(mean_d(zk),var_d(zk))
%              w ~ N(0,sigmaw)
%
%   
%   x = [s, ds, th, dth]'   carriage position and pole angle (and derivatives)
%   u = [F]'                force on the carriage and torque on the pole joint
%   
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% TODO:
%       -[ ] Implement the case for multiple inputs
%--------------------------------------------------------------------------

    properties
        Mc      % mass of the carriage
        Mp      % mass of the pole
        b       % friction coefficient between the carriage and the floor
        I       % inertia matrix of the pole CG
        l       % pole length
        
        d       % [E[d(z)] , Var[d(z)]] = d(z): disturbace model
        w       % [E[w] , Var[w]] = w: measurement noise
        
        g = 9.8;
    end
    
    properties(SetAccess=private)
        Bd = [0;            % xk+1 = fd(xk,uk) + Bd*d(zk)
              0;
              1;
              0]
        Bz = [0 0 1 0       % z = Bz*x
              0 0 0 1];     
        n = 4               % number of outputs x(t)
        m = 1               % number of inputs u(t)
        nd                  % output dimension of du(z)
        nz                  % dimension of z(t)
    end
    
    
    methods
        
        function obj = InvertedPendulumModel (Mc, Mp, b, I, l, d, sigmaw)
        %------------------------------------------------------------------
        %   object constructor
        %------------------------------------------------------------------
            obj.Mc = Mc;
            obj.Mp = Mp;
            obj.b = b;
            obj.I = I;
            obj.l = l;
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
        %   dimension of z(t)
        %------------------------------------------------------------------
            nz = size(obj.Bz,1);
        end
        
        function [xdot, grad_xdot] = f (obj, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics of the inverted pendulum (including
        %   disturbance):
        %
        %       p(xdot | x,u) = N(xdot, grad_xdot'*var_xdot*grad_xdot)
        %
        %   (Mc+Mp)*dds + b*ds + Mp*l/2*ddth*cos(th) - Mp*l/2*dth^2*sin(th) = F
        %   (I+Mp*(l/2)^2)*ddth + Mp*g*l/2*sin(th) + Mp*l*dds*cos(th) = T
        %   
        %   x = [s, ds, th, dth]'
        %   u = [F]'
        %
        %   (HOW TO GENERATE THESE EQUATIONS:)
        %
        % syms g Mc Mp b I l F T s ds dds  th dth ddth real
        % fzero = [(Mc+Mp)*dds + b*ds + Mp*l/2*ddth*cos(th) - Mp*l/2*dth^2*sin(th) - F ;
        %        (I+Mp*(l/2)^2)*ddth + Mp*g*l/2*sin(th) + Mp*l*dds*cos(th) - T  ];
        % sol = solve(fzero,[dds,ddth])
        % dds = simplify(sol.dds)
        % ddth = simplify(sol.ddth)
        % 
        % u = [F,T]';
        % x = [s, ds, th, dth]'
        % xdot = [ds, dds, dth, ddth]';
        % params = [Mc Mp I g l b ]';
        % 
        % grad_xdot = simplify(jacobian(xdot,x)');
        % matlabFunction( xdot, grad_xdot, 'Vars', {x;u;params} ,'File', 'invertedPendulumModel_f' )
        %------------------------------------------------------------------
            
            params = [obj.Mc obj.Mp obj.I obj.g obj.l obj.b]';
            u(2,1) = 0;  % set Torque to zero
            [xdot, grad_xdot] = invertedPendulumModel_f( x, u, params );
        end
        
        
        function [xkp1, grad_xkp1] = fd (obj, xk, uk, dt)
        %------------------------------------------------------------------
        %   Discrete time dynamics of the inverted pendulum (including
        %   disturbance)
        %------------------------------------------------------------------
            % calculate continous time dynamics
            [xdot, grad_xdot] = obj.f(xk,uk);
            
            % discretize
            xkp1      = xk + dt * xdot;
            grad_xkp1 = eye(obj.n) + dt * grad_xdot;
        end
        
        
        function [mu_xkp1,var_xkp1] = xkp1 (obj, mu_xk, var_xk, uk, dt)
        %------------------------------------------------------------------
        %   State prediction (motion model) of inverted pendulum using
        %   Extended Kalman Filter approach
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
            var_x_d_w = blkdiag(var_xk, var_d, var_w);
            
            % predict mean and variance (Extended Kalman Filter)
            mu_xkp1  = fd  + obj.Bd * ( mu_d + mu_w );
            var_xkp1 = grad_xkp1' * var_x_d_w * grad_xkp1;
        end
        
        
        function [A,B] = linearize (obj)
        %------------------------------------------------------------------
        % Return continuous time linearized model parameters A,B
        %       xdot = A*x + B*u
        %------------------------------------------------------------------
            Mc=obj.Mc; Mp=obj.Mp; b=obj.b; I=obj.I; l=obj.l; g=obj.g;
            p = I*(Mc+Mp)+Mc*Mp*l^2;
            A = [0      1              0           0;
                 0 -(I+Mp*l^2)*b/p  (Mp^2*g*l^2)/p   0;
                 0      0              0           1;
                 0 -(Mp*l*b)/p       Mp*g*l*(Mc+Mp)/p  0];
            B = [     0;
                 (I+Mp*l^2)/p;
                      0;
                    Mp*l/p];
        end
        
    end
end
