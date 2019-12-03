classdef invertedPendulum
%------------------------------------------------------------------
%   dot_x(t) = f_true(x(t),u(t)) + Bd*(d_true(t) + w(t)),    wk=N(0,sigmaw^2)
%   
%   x = x = [x, dx, th, dth]'
%   u = [F, T]'
%------------------------------------------------------------------
    
    properties
        M = 0.5;
        m = 0.2;
        b = 0.1;
        I = 0.006;
        g = 9.8;
        l = 0.6;
        sigmaw = 0.01;
        
        Bd
    end
    
    methods
        
        function obj = invertedPendulum()
        end
        
        
        function d = d_true(obj,x,u)
            d = [20 -20]*[mvnpdf([x,u],[2,2],diag([2,20])) mvnpdf([x,u],[-2,-2],diag([2,20]))]';
        end
        
        
        function xdot = f_true(obj,x,u)
        %------------------------------------------------------------------
        %   (M+m)*ddx + b*dx + m*l/2*ddth*cos(th) - m*l/2*dth^2*sin(th) = F
        %   (I+m*(l/2)^2)*ddth + m*g*l/2*sin(th) + m*l*ddx*cos(th) = T
        %   
        %   x = [x, dx, th, dth]'
        %   u = [F, T]'
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
            x  = x(1);
            dx = x(2);
            th = x(3);
            dth= x(4);
            
            F = u(1);
            T = u(2);
            
            % simplify code notation
            M=obj.M; m=obj.m; I=obj.I; g=obj.g; l=obj.l; b=obj.b;
            % equations of motion
            ddx  = (8*F*I - 8*I*b*dx + 2*F*l^2*m - 2*b*dx*l^2*m - 4*T*l*m*cos(th) + dth^2*l^3*m^2*sin(th) + 2*g*l^2*m^2*cos(th)*sin(th) + 4*I*dth^2*l*m*sin(th))/(2*((1 - 2*cos(th)^2)*l^2*m^2 + M*l^2*m + 4*I*m + 4*I*M));
            ddth = -(2*(2*F*l*m*cos(th) - 2*T*m - 2*M*T + g*l*m^2*sin(th) + dth^2*l^2*m^2*cos(th)*sin(th) + M*g*l*m*sin(th) - 2*b*dx*l*m*cos(th)))/((1 - 2*cos(th)^2)*l^2*m^2 + M*l^2*m + 4*I*m + 4*I*M);
            
            xdot = [dx, ddx, dth, ddth]';
        end
        
        
        function xkp1 = fd_true(obj,xk,uk,dt,inclNoise)
        %------------------------------------------------------------------
        %   Discrete time dynamics of the inverted pendulum
        %------------------------------------------------------------------
            xkp1 = xk + dt * ( obj.f_true(obj,xk,uk) + obj.Bd*obj.d_true(x,u));
            
            if inclNoise
                dicreteNoiseVariance = dt * obj.sigmaw^2;
                wk = sqrt(dicreteNoiseVariance) * randn(size(obj.Bd,2));
                xkp1 = xkp1 + obj.Bd * wk;
            end
            
        end
       
        
    end
end






