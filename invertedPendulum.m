classdef invertedPendulum
%------------------------------------------------------------------
%   dot_x(t) = f_true(x(t),u(t)) + Bd*(du(z(t)) + w(t)),    
%
%       where: z(t)=Bz*x(t) and wk=N(0,sigmaw^2)
%   
%   x = [x, dx, th, dth]'   carriage position and pole angle (and derivatives)
%   u = [F, T]'             force on the carriage and torque on the pole joint
%   
%   du = unmodeled disturbances
%------------------------------------------------------------------
    
    properties
        Mc      % mass of the carriage
        Mp      % mass of the pole
        b       % friction coefficient between the carriage and the floor
        I       % inertia matrix of the pole CG
        l       % pole length
        sigmaw  % du(t) std deviation
        
        g = 9.8;
    end
    
    properties(SetAccess=private)
        Bd = [0;            % unmodeled dynamics only affect velocity components
              0;
              1; 
              0]
        Bz = [0 1 0 0       % z = Bz*x
              0 0 0 1];     
        n = 2       % number of outputs x(t)
        m = 2       % number of inputs u(t)
        nd          % output dimension of du(z)
        nz          % dimension of z(t)
    end
    
    
    methods
        
        function obj = invertedPendulum (Mc, Mp, b, I, l, sigmaw)
            obj.Mc = Mc;
            obj.Mp = Mp;
            obj.b = b;
            obj.I = I;
            obj.l = l;
            obj.sigmaw = sigmaw;
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
        
        
        function d = du(~,z,u)
            %------------------------------------------------------------------
            %   Unmodeled disturbances
            % args
            %   d: nz->nd
            %------------------------------------------------------------------
            d = [20 -20]*[mvnpdf([z,u],[2,2],diag([2,20])) mvnpdf([z,u],[-2,-2],diag([2,20]))]';
        end
        
        
        function xdot = f (obj,x,u)
        %------------------------------------------------------------------
        %   Continuous time dynamics of the inverted pendulum (including
        %   disturbance)
        %
        %   (Mc+Mp)*ddx + b*dx + Mp*l/2*ddth*cos(th) - Mp*l/2*dth^2*sin(th) = F
        %   (I+Mp*(l/2)^2)*ddth + Mp*g*l/2*sin(th) + Mp*l*ddx*cos(th) = T
        %   
        %   x = [th, dth]'
        %   u = [F, T]'
        %   d = [x, dx]'
        %
        %   (HOW TO GENERATE THIS EQUATIONS:)
        %
        %   syms g Mc Mp b I l F T x dx ddx  th dth ddth
        %   fzero = [(Mc+Mp)*ddx + b*dx + Mp*l/2*ddth*cos(th) - Mp*l/2*dth^2*sin(th) - F ;
        %            (I+Mp*(l/2)^2)*ddth + Mp*g*l/2*sin(th) + Mp*l*ddx*cos(th) - T  ];
        %   sol = solve(fzero,[ddx,ddth])
        %   ddx = simplify(sol.ddx)
        %   ddth = simplify(sol.ddth)
        %   matlabFunction([ddx;ddth],'file','invertedPendulum_f_true')
        %------------------------------------------------------------------
            x  = x(1);
            dx = x(2);
            th = x(3);
            dth= x(4);
            
            F = u(1);
            T = u(2);
            
            % simplify code notation
            Mc=obj.Mc; Mp=obj.Mp; I=obj.I; g=obj.g; l=obj.l; b=obj.b;
            
            % equations of motion
            ddx  = (8*F*I - 8*I*b*dx + 2*F*l^2*m - 2*b*dx*l^2*m - 4*T*l*m*cos(th) + dth^2*l^3*m^2*sin(th) + 2*g*l^2*m^2*cos(th)*sin(th) + 4*I*dth^2*l*m*sin(th))/(2*((1 - 2*cos(th)^2)*l^2*m^2 + M*l^2*m + 4*I*m + 4*I*M));
            ddth = -(2*(2*F*l*Mp*cos(th) - 2*T*Mp - 2*Mc*T + g*l*Mp^2*sin(th) + dth^2*l^2*Mp^2*cos(th)*sin(th) + Mc*g*l*Mp*sin(th) - 2*b*dx*l*Mp*cos(th)))/((1 - 2*cos(th)^2)*l^2*Mp^2 + Mc*l^2*Mp + 4*I*Mp + 4*I*Mc);
            
            xdot = [dx, ddx, dth, ddth]';
        end
        
        
        function xkp1 = fd (obj,xk,uk,dt,inclNoise)
        %------------------------------------------------------------------
        %   Discrete time dynamics of the inverted pendulum (including
        %   disturbance)
        %------------------------------------------------------------------
            zk = obj.Bz * xk; 
            xkp1 = xk + dt * ( obj.f(obj,xk,uk) + obj.Bd*obj.du(zk,u));
            
            if inclNoise
                dicreteNoiseVariance = dt * obj.sigmaw^2;
                wk = sqrt(dicreteNoiseVariance) * randn(obj.md);
                xkp1 = xkp1 + obj.Bd * wk;
            end    
        end
       
        
    end
end






