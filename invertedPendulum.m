classdef invertedPendulum
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
        
        d       % [mean_dz,var_dz] = d(z): disturbace model
        
        g = 9.8;
    end
    
    properties(SetAccess=private)
        Bd = [0;            % how unmodeled dynamics affect states
              0;
              1;
              0]
        Bz = [0 0 1 0       % z = Bz*x
              0 0 0 1];     
        n = 4       % number of outputs x(t)
        m = 1       % number of inputs u(t)
        nd          % output dimension of du(z)
        nz          % dimension of z(t)
    end
    
    
    methods
        
        function obj = invertedPendulum (Mc, Mp, b, I, l, d)
        %------------------------------------------------------------------
        %   object constructor
        %------------------------------------------------------------------
            obj.Mc = Mc;
            obj.Mp = Mp;
            obj.b = b;
            obj.I = I;
            obj.l = l;
            obj.d = d;
        end
        
        
        function nd = get.nd(obj)
        %------------------------------------------------------------------
        %   output dimension of du(z)
        %------------------------------------------------------------------
            nd = size(obj.Bd,2);
        end
        
        function nz = get.nz(obj)
        %------------------------------------------------------------------
        %   dimension of z(t)
        %------------------------------------------------------------------
            nz = size(obj.Bz,1);
        end
        
        function [mu_xdot, var_xdot] = f (obj,x,u)
        %------------------------------------------------------------------
        %   Continuous time dynamics of the inverted pendulum (including
        %   disturbance):
        %
        %       p(xdot | x,u) = N(mu_xdot,var_xdot)
        %
        %   (Mc+Mp)*dds + b*ds + Mp*l/2*ddth*cos(th) - Mp*l/2*dth^2*sin(th) = F
        %   (I+Mp*(l/2)^2)*ddth + Mp*g*l/2*sin(th) + Mp*l*dds*cos(th) = T
        %   
        %   x = [s, ds, th, dth]'
        %   u = [F]'
        %
        %   (HOW TO GENERATE THIS EQUATIONS:)
        %
        %   syms g Mc Mp b I l F T s ds dds  th dth ddth
        %   fzero = [(Mc+Mp)*ddx + b*ds + Mp*l/2*ddth*cos(th) - Mp*l/2*dth^2*sin(th) - F ;
        %            (I+Mp*(l/2)^2)*ddth + Mp*g*l/2*sin(th) + Mp*l*dds*cos(th) - T  ];
        %   sol = solve(fzero,[dds,ddth])
        %   dds = simplify(sol.dds)
        %   ddth = simplify(sol.ddth)
        %   matlabFunction([dds;ddth],'file','invertedPendulum_f_true')
        %------------------------------------------------------------------
            s  = x(1);
            ds = x(2);
            th = x(3);
            dth= x(4);
            
            F = u(1);
            T = 0;
            
            % simplify code notation
            Mc=obj.Mc; Mp=obj.Mp; I=obj.I; g=obj.g; l=obj.l; b=obj.b;
            
            % equations of motion
            th = th + pi;
            dds  = (8*F*I - 8*I*b*ds + 2*F*l^2*Mp - 2*b*ds*l^2*Mp - 4*T*l*Mp*cos(th) + dth^2*l^3*Mp^2*sin(th) + 2*g*l^2*Mp^2*cos(th)*sin(th) + 4*I*dth^2*l*Mp*sin(th))/(2*((1 - 2*cos(th)^2)*l^2*Mp^2 + Mc*l^2*Mp + 4*I*Mp + 4*I*Mc));
            ddth = -(2*(2*F*l*Mp*cos(th) - 2*T*Mp - 2*Mc*T + g*l*Mp^2*sin(th) + dth^2*l^2*Mp^2*cos(th)*sin(th) + Mc*g*l*Mp*sin(th) - 2*b*ds*l*Mp*cos(th)))/((1 - 2*cos(th)^2)*l^2*Mp^2 + Mc*l^2*Mp + 4*I*Mp + 4*I*Mc);
            
            % calculate xdot mean and covariance
            mu_xdot  = [ds, dds, dth, ddth]';
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






