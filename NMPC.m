%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -
%------------------------------------------------------------------

classdef NMPC < handle
    %------------------------------------------------------------------
    % Nonlinear MPC class
    % 
    %   solve:
    %
    %   MIN     { SUM_{k=0:N-1} fo(tk,xk,uk,r(t)) } + fend(tN,xN,r(tN))
    %
    %   s.t.    xk+1 = E[f(xk,uk)]
    %           h(xk,uk) == 0
    %           g(xk,uk) <= 0
    %
    %   where the motion model evaluates   [E[xk+1],Var[xk+1]] = f(xk,uk)
    %
    %   for x0, u0,...,uN-1, e0,...,eN
    %
    %   where xk: state variables
    %         zk: selected state variables zk=Bd'*xk
    %         r(tk): trajectory
    %         tk: current time
    %------------------------------------------------------------------
    
    properties
        % Optimizer settings
        maxiter = 200       % max solver iterations
        tol     = 1e-6      % optimizer tolerance
        N       = 30        % prediction horizon
        
        % Define optimization problem
        f       % @fun nonlinear dynamics:  [E[xk+1],Var[xk+1]] = f(xk,uk)
        h       % equality constraint function 
        g       % inequality constraint function
        
        fo      % @fun nonlinear cost function
        fend    % @fend nonlinear cost function for the final state
        
        % Optimization dimension
        n   % dim(x) state space dimension
        m   % dim(u) input diension
        dt  % time step size
        nh  % number of additional eq. constraints for every time step
        ng  % number of additional ineq. constraints for every time step
    end
    
    properties(Access=private)
        % save last optimal results computed, in order to use as initial guess
        vars_opt_old = []
    end
    
    
    methods
        
        function obj = NMPC (f, h, g, n, m, fo, fend, N, dt)
        %------------------------------------------------------------------
        % MPC constructor
        %
        %   f: motion model that evaluates  [E[xk+1],Var[xk+1]] = f(xk)
        %------------------------------------------------------------------
           % constraints
           obj.f  = f;
           obj.h = h;
           obj.g = g;
           % get size of additional constraints
           obj.nh = length(h(zeros(n,1),zeros(m,1)));
           obj.ng = length(h(zeros(n,1),zeros(m,1)));
           
           % variable dimensions
           obj.n = n;
           obj.m = m;
           % cost functions
           obj.fo = fo;
           obj.fend = fend;
           % optimizer parameters
           obj.N = N;
           obj.dt = dt;
        end
        
        
        function numvars = optSize(obj)
        %------------------------------------------------------------------
        % How many variables we need to optimize?
        %
        %   vars_opt = [u0,...,uN-1]
        %------------------------------------------------------------------
            numvars = obj.n + obj.N*obj.m;
        end
        
        
        function u0_opt = optimize(obj, x0, t0, r)
        %------------------------------------------------------------------
        % Calculate first uptimal control input
        %------------------------------------------------------------------
            
            %-------- Set initial guess for optimization variables  -------
            % initialize optimization variables initial guesses
            if isempty(obj.vars_opt_old)  % if this is the first optimization
                uguess = zeros(obj.m * obj.N,1);      % [ u0,...,uN-1]
            else
                [~,uguess] = obj.splitvariables(obj.vars_opt_old);
                uguess = uguess(:); % <m,N> to <m+N,1>
            end
            varsguess = [x0; uguess];
            %--------------------------------------------------------------
            
            assert( numel(varsguess) == obj.optSize(), ...
                'There is something wrong with the code. Number of optimization variables does not match!' );
                
            % We currently dont need this parameters in fmincon
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            lb = [];
            ub = [];
            
            % define cost and constr. function, depending only on the
            % variables. This is a prerequisite for the function fmincon
            costfun = @(vars) obj.costfun(vars, t0, r);
            nonlcon = @(vars) obj.nonlcon(vars, t0, x0);
            
            % define optimizer settings
            options = optimoptions('fmincon',...
                                   'Display','iter',...
                                   'Algorithm','sqp',...
                                   'ConstraintTolerance',obj.tol,...
                                   'MaxIterations',obj.maxiter);
            
            % solve optimization problem                   
            vars_opt = fmincon(costfun,varsguess,A,b,Aeq,beq,lb,ub,nonlcon,options);

            % store current optimization results to use as initial guess
            % for future optimizations
            obj.vars_opt_old = vars_opt;
            
            % split variables since vars_opt = [x_opt; u_opt; e_opt]
            [~, u_opt] = splitvariables(obj, vars_opt);
            u0_opt = u_opt(:,1);
        end
        
    
        function [x0, uvec] = splitvariables(obj, vars)
        %------------------------------------------------------------------
        % args:
        %   vars: <optSize,1> optimization variables
        % out:
        %   x0: <n,1>
        %   uvec: <m,N>
        %------------------------------------------------------------------
            % split variables
            x0   = vars(1:obj.n);
            uvec = vars( (1:obj.N*obj.m)  + length(x0) );
            % reshape the column vector to <m,N>
            uvec = reshape(uvec, obj.m, obj.N);
        end
        
        
        function [mu_xk,var_xk] = calculateStateSequence(obj, mu_x0, var_x0, uk)
        %------------------------------------------------------------------
        % Propagate mean and covariance of state sequence, given control
        % input sequence.
        %------------------------------------------------------------------
            mu_xk  = zeros(obj.n,obj.N+1);
            var_xk = zeros(obj.n,obj.n,obj.N+1);
            
            mu_xk(:,1) = mu_x0;
            var_xk(:,:,1) = var_x0;
            
            for iN=1:obj.N      % [x1,...,xN]
                [mu_xk(:,iN+1),var_xk(:,:,iN+1)] = obj.f(mu_xk(:,iN),var_xk(:,:,iN),uk(:,iN));
            end
        end

        
        function cost = costfun(obj, vars, t0, r)
        %------------------------------------------------------------------
        % Evaluate cost function for the whole horizon, given variables
        %------------------------------------------------------------------
            % split variables
            [mu_x0, uk] = obj.splitvariables(vars);
            var_x0 = zeros(obj.n);
            
            % calculate state sequence for given control input sequence and x0
            [mu_xk,var_xk] = obj.calculateStateSequence(mu_x0, var_x0, uk);
                        
            cost = 0;
            t = t0;
            for iN=1:obj.N      % i=0:N-1
                % add cost
                cost = cost + obj.fo(t, mu_xk(:,iN), var_xk(:,:,iN), uk(:,iN), r);

                % update current time
                t = t + iN * obj.dt;
            end
            % final cost
            cost = cost + obj.fend(t, mu_xk(:,end), var_xk(:,:,end), r);
        end
        

        function [cineq,ceq] = nonlcon(obj, vars, t0, x0)
        %------------------------------------------------------------------
        % Evaluate nonlinear equality and inequality constraints
        % args
        %   cineq = g(x,u)
        %   ceq   = h(x,u)
        %------------------------------------------------------------------ 

            % init vectors to speedup calculations
            ceq_dyn = zeros(obj.n,  1);
            ceq_h   = zeros(obj.nh, obj.N);
            cineq_g = zeros(obj.ng, obj.N);
        
            % split variables
            [mu_x0, uk] = obj.splitvariables(vars);
            var_x0 = zeros(obj.n);
            
            % calculate state sequence for given control input sequence and x0
            [mu_xk,var_xk] = obj.calculateStateSequence(mu_x0, var_x0, uk);
            
            % set initial state constraint: x0 - x(0) = 0
            ceq_dyn(:,1) = x0 - mu_xk(:,1);
            
            t = t0;
            for iN=1:obj.N
                
                % append provided equality constraints(h==0)
                ceq_h(:,iN) = obj.h(mu_xk(:,iN),uk(:,iN));
                
                % provided inequality constraints (g<=0)
                cineq_g(:,iN) = obj.g(mu_xk(:,iN),uk(:,iN));

                t = t + iN * obj.dt;
            end
            ceq   = [ceq_dyn(:); ceq_h(:)];
            cineq = cineq_g(:);
        end
        
    end
end


