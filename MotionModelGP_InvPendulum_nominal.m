%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -
%------------------------------------------------------------------


classdef MotionModelGP_InvPendulum_nominal < MotionModelGP
%--------------------------------------------------------------------------
%   xk+1 = fd(xk,uk) + Bd * ( d(zk) + w ),    
%
%       where: zk = [Bz_x*xk ; Bz_u*uk],
%              d ~ N(mean_d(zk),var_d(zk))
%              w ~ N(0,sigmaw)
%
%   
%   x = [s, ds, th, dth]'   carriage position and pole angle (and derivatives)
%   u = [F]'                force on the carriage and torque on the pole joint
%   
%--------------------------------------------------------------------------

    properties
        Mc      % mass of the carriage
        Mp      % mass of the pole
        b       % friction coefficient between the carriage and the floor
        I       % inertia matrix of the pole CG
        l       % pole length
        g = 9.8
    end
    
    properties(Constant)
        % keep in mind the dimensions:  xk+1 = fd(xk,uk) + Bd*(d(z)+w)),
        % where z = [Bz_x*x;Bz_u*u] 
        Bz_x = [0 0 1 0
                0 0 0 1] 
        Bz_u = []; 
        Bd = [0;            % xk+1 = fd(xk,uk) + Bd*d(zk)
              0;
              1;
              0]
            
        n  = 4   % number of outputs x(t)
        m  = 1   % number of inputs u(t)
        nz = 2   % dimension of z(t)
        nd = 1   % output dimension of d(z)
    end
    
    
    methods
        
        function obj = MotionModelGP_InvPendulum_nominal (Mc, Mp, b, I, l, d, sigmaw)
        %------------------------------------------------------------------
        %   object constructor
        %------------------------------------------------------------------
            % call superclass constructor
            obj = obj @ MotionModelGP(d,sigmaw);
            % store parameters
            obj.Mc = Mc;
            obj.Mp = Mp;
            obj.b = b;
            obj.I = I;
            obj.l = l;
            
            % add folder CODEGEN to path. Here there will be some functions
            % generated with the method generate_invertedPendulum_functions()
            addpath(fullfile(pwd,'CODEGEN'))
        end
        
        function xdot = f (obj, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics.
        %   out:
        %       xdot: <n,1> time derivative of x given x and u
        %------------------------------------------------------------------
            params = [obj.Mc obj.Mp obj.I obj.g obj.l obj.b]';
            xdot = invertedPendulum_f(x, u, params );
        end
        
        function gradx = gradx_f(obj, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics.
        %   out:
        %       gradx: <n,n> gradient of xdot w.r.t. x
        %------------------------------------------------------------------
            params = [obj.Mc obj.Mp obj.I obj.g obj.l obj.b]';
            gradx = invertedPendulum_gradx_f(x, u, params );
        end
        
        function gradu = gradu_f(obj, x, u)
        %------------------------------------------------------------------
        %   Continuous time dynamics.
        %   out:
        %       gradu: <m,n> gradient of xdot w.r.t. u
        %------------------------------------------------------------------
            params = [obj.Mc obj.Mp obj.I obj.g obj.l obj.b]';
            gradu = invertedPendulum_gradu_f(x, u, params );
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
    
    methods(Static)
        function generate_invertedPendulum_functions()
        %------------------------------------------------------------------
        %   Generate continuous time dynamics equations of the inverted 
        %   pendulum:. This function generates three functions:
        %       xdot = f(x,u)       - dynamics
        %       gradx_xdot(x,u)     - gradient of xdot w.r.t. x
        %       gradu_xdot(x,u)     - gradient of xdot w.r.t. u
        %
        %   (Mc+Mp)*dds + b*ds + Mp*l/2*ddth*cos(th) - Mp*l/2*dth^2*sin(th) = F
        %   (I+Mp*(l/2)^2)*ddth + Mp*g*l/2*sin(th) + Mp*l*dds*cos(th) = T
        % 
        %   x = [s, ds, th, dth]'
        %   u = [F]'
        %
        %
        %   Example how to run this function:
        %       ip = MotionModelGP_InvertedPendulum(5, 2, 0.1, 0.6, 3, @(z)deal(0,0), 0);
        %       ip.generate_invertedPendulum_functions();
        %------------------------------------------------------------------
            syms g Mc Mp b I l F T s ds dds  th dth ddth real
            T = 0;  % we are not using this input for now
            fzero = [(Mc+Mp)*dds + b*ds + Mp*l/2*ddth*cos(th) - Mp*l/2*dth^2*sin(th) - F ;
                   (I+Mp*(l/2)^2)*ddth + Mp*g*l/2*sin(th) + Mp*l*dds*cos(th) - T  ];
            sol = solve(fzero,[dds,ddth]);
            dds = simplify(sol.dds);
            ddth = simplify(sol.ddth);

            u = F;
            x = [s, ds, th, dth]';
            xdot = [ds, dds, dth, ddth]';
            params = [Mc Mp I g l b ]';
            
            
            folder = fullfile(pwd,'CODEGEN');
            if ~exist(folder,'dir')
                mkdir(folder); 
            end
            addpath(folder)
            
            
            matlabFunction( xdot, 'Vars', {x;u;params} ,'File', fullfile('CODEGEN','invertedPendulum_f') );

            gradx_f = simplify(jacobian(xdot,x)');
            matlabFunction( gradx_f, 'Vars', {x;u;params} ,'File', fullfile('CODEGEN','invertedPendulum_gradx_f') );

            gradu_f = simplify(jacobian(xdot,u)');
            matlabFunction( gradu_f, 'Vars', {x;u;params} ,'File', fullfile('CODEGEN','invertedPendulum_gradu_f') );
        
            disp('FINISHED! functions invertedPendulum_f, invertedPendulum_gradx_f and invertedPendulum_gradu_f generated!!')

        end
    end
    
end
