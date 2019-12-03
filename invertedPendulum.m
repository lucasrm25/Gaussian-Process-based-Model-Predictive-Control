classdef invertedPendulum
%------------------------------------------------------------------
%   dot_x(t) = f_true(x(t),u(t),d(t)) + Bd*(du(t) + w(t)),    wk=N(0,sigmaw^2)
%   
%   x = [th, dth]'   angle and angle derivative of pole angle
%   u = [F, T]'      force on the carriage and torque on the pole joint
%   d = [x, dx]'     carriage position and velocity
%
%   du = unmodeled disturbances
%------------------------------------------------------------------
    
    properties
        M       % mass of the carriage
        m       % mass of the pole
        b       % friction coefficient between the carriage and the floor
        I       % inertia matrix of the pole CG
        l       % pole length
        sigmaw  % du(t) std deviation
        
        g = 9.8;
        
        Bd
    end
    
    methods
        
        function obj = invertedPendulum (M, m, b, I, l, Bd, sigmaw)
            obj.M = M;
            obj.m = m;
            obj.b = b;
            obj.I = I;
            obj.l = l;
            obj.Bd = Bd;
            obj.sigmaw = sigmaw;
        end
        
        
        function d = du(~,x,u)
            %------------------------------------------------------------------
            %   Unmodeled disturbances
            %------------------------------------------------------------------
            d = [20 -20]*[mvnpdf([x,u],[2,2],diag([2,20])) mvnpdf([x,u],[-2,-2],diag([2,20]))]';
        end
        
        
        function xdot = f (obj,x,u,d)
        %------------------------------------------------------------------
        %   Continuous time dynamics of the inverted pendulum (including
        %   disturbance)
        %
        %   (M+m)*ddx + b*dx + m*l/2*ddth*cos(th) - m*l/2*dth^2*sin(th) = F
        %   (I+m*(l/2)^2)*ddth + m*g*l/2*sin(th) + m*l*ddx*cos(th) = T
        %   
        %   x = [th, dth]'
        %   u = [F, T]'
        %   d = [x, dx]'
        %
        %   (HOW TO GENERATE THIS EQUATIONS:)
        %
        %   syms g m M b I l F T x dx ddx  th dth ddth
        %   fzero = [(M+m)*ddx + b*dx + m*l/2*ddth*cos(th) - m*l/2*dth^2*sin(th) - F ;
        %          (I+m*(l/2)^2)*ddth + m*g*l/2*sin(th) + m*l*ddx*cos(th) - T  ];
        %   sol = solve(fzero,[ddx,ddth])
        %   ddx = simplify(sol.ddx)
        %   ddth = simplify(sol.ddth)
        %   matlabFunction([ddx;ddth],'file','invertedPendulum_f_true')
        %------------------------------------------------------------------
            th = x(1);
            dth= x(2);
            
            F = u(1);
            T = u(2);
            
            x  = d(1);
            dx = d(2);
            
            % simplify code notation
            M=obj.M; m=obj.m; I=obj.I; g=obj.g; l=obj.l; b=obj.b;
            
            % equations of motion
            % ddx  = (8*F*I - 8*I*b*dx + 2*F*l^2*m - 2*b*dx*l^2*m - 4*T*l*m*cos(th) + dth^2*l^3*m^2*sin(th) + 2*g*l^2*m^2*cos(th)*sin(th) + 4*I*dth^2*l*m*sin(th))/(2*((1 - 2*cos(th)^2)*l^2*m^2 + M*l^2*m + 4*I*m + 4*I*M));
            ddth = -(2*(2*F*l*m*cos(th) - 2*T*m - 2*M*T + g*l*m^2*sin(th) + dth^2*l^2*m^2*cos(th)*sin(th) + M*g*l*m*sin(th) - 2*b*dx*l*m*cos(th)))/((1 - 2*cos(th)^2)*l^2*m^2 + M*l^2*m + 4*I*m + 4*I*M);
            
            xdot = [dx, ddx, dth, ddth]';
        end
        
        
        function xkp1 = fd (obj,xk,uk,dk,dt,inclNoise)
        %------------------------------------------------------------------
        %   Discrete time dynamics of the inverted pendulum (including
        %   disturbance)
        %------------------------------------------------------------------
            xkp1 = xk + dt * ( obj.f(obj,xk,uk,dk) + obj.Bd*obj.du(x,u));
            
            if inclNoise
                dicreteNoiseVariance = dt * obj.sigmaw^2;
                wk = sqrt(dicreteNoiseVariance) * randn(size(obj.Bd,2));
                xkp1 = xkp1 + obj.Bd * wk;
            end    
        end
       
        
    end
end






