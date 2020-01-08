classdef ODE < handle
    
    properties(Constant)
        RK4 = [ 0 0 0 0 0;
                .5 .5 0 0 0;
                .5 0 .5 0 0;
                1 0 0 1 0;
                0 1/6 1/3 1/3 1/6]
        RK2 = [0 0 0; 
               .5 .5 0; 
               0 0 1]
        RK1 = [0 0; 
               0 1]
    end
    
    methods
        function obj = ODE()
        end
    end
    
    methods(Static)
        
        function [t,x] = ode1(f, x0, t_end, dt)
            [t,x] = ODE.RK(ODE.RK1, f, x0, t_end, dt);
        end
        
        function [t,x] = ode2(f, x0, t_end, dt)
            [t,x] = ODE.RK(ODE.RK2, f, x0, t_end, dt);
        end
        
        function [t,x] = ode4(f, x0, t_end, dt)
            [t,x] = ODE.RK(ODE.RK4, f, x0, t_end, dt);
        end
        
        function [t,x] = RK(ButcherT, f, x0, t_end, dt)
            %------------------------------------------------------------------
            % ButcherT: <n,n> Butcher tableau of order n-1
            % f: <@(t,x)> ode function, e.g. f = @(t,x) lambda*x;
            %------------------------------------------------------------------
            dim = size(x0,1);
            ord = size(ButcherT,1)-1;
            N = ceil(t_end/dt);
            t = zeros(1,N);
            x = zeros(dim,N);
            x(:,1) = x0;    
            for it = 2:N+1
                t(it) = t(it-1) + dt;
                K = zeros(ord,dim);
                for ki=1:ord
                    K(ki,:) = f( t(it-1)+ButcherT(ki,1)*dt , x(:,it-1)+ dt*(ButcherT(ki,2:end)*K)' );
                end      
                x(:,it) = x(:,it-1) + dt * (ButcherT(end,2:end)*K)';
            end
            t = t(2:end);
            x = x(:,2:end);
        end
    end
end

