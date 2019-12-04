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
    %   for x0,...,xN, u0,...,uN-1, e0,...,eN
    %
    %   where xk: state variables
    %         zk: selected state variables zk=Bd'*xk
    %         ek: extra variables
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
        %   vars_opt = [x0,...,xN, u0,...,uN-1, e0,...,eN]
        %------------------------------------------------------------------
            numvars = (obj.N+1)*(obj.n) + (obj.N)*(obj.m);
        end
        
        
        function u0_opt = optimize(obj, x0, t0, r)
        %------------------------------------------------------------------
        % Calculate first uptimal control input
        %------------------------------------------------------------------
            
            %-------- Set initial guess for optimization variables  -------
            % initialize optimization variables initial guesses
            if isempty(obj.vars_opt_old)
                % if this is the first optimization
                xguess =  repmat(x0,obj.N+1,1);  % [ x0,...,xN-1,xN ]
                uguess = zeros(obj.m * obj.N,1);           % [ u0,...,uN-1    ]
                % vector of initial guesses
                varsguess = [xguess; uguess];
            else
                varsguess = obj.vars_opt_old;
                varsguess(1:obj.n) = x0;
            end
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
        
    
        function [xvec, uvec] = splitvariables(obj, vars)
            % split variables
            xvec = vars(1:(obj.N+1)*obj.n);
            uvec = vars( (1:obj.N*obj.m)  + length(xvec) );
            % reshape the column vector to <n,N+1> and <m,N>
            xvec = reshape(xvec, obj.n, obj.N+1);
            uvec = reshape(uvec, obj.m, obj.N);
        end
        

        function cost = costfun(obj, vars, t0, r)
        %------------------------------------------------------------------
        % Evaluate cost function for the whole horizon, given variables
        %------------------------------------------------------------------
            % split variables
            [xk, uk] = obj.splitvariables(vars);

            cost = 0;
            t = t0;
            for iN=1:obj.N      % i=0:N-1
                % add cost
                cost = cost + obj.fo(t,xk(:,iN),uk(:,iN),r);
                
                % update current time
                t = t + iN * obj.dt;
            end
            % final cost
            cost = cost + obj.fend(t,xk(:,end-obj.n+1),r);
        end
        

        function [cineq,ceq] = nonlcon(obj,vars, t0, x0)
        %------------------------------------------------------------------
        % Evaluate nonlinear equality and inequality constraints
        % args
        %   cineq = g(x,u)
        %   ceq   = h(x,u)
        %------------------------------------------------------------------ 

            % init vectors to speedup calculations
            ceq_dyn = zeros(obj.n,  obj.N+1);
            ceq_h   = zeros(obj.nh, obj.N);
            cineq_g = zeros(obj.ng, obj.N);
        
            % split variables
            [xk, uk] = obj.splitvariables(vars);
            
            % set initial state constraint: x0 - x(0) = 0
            ceq_dyn(:,1) = x0 - xk(:,1);
            
            t = t0;            
            for iN=1:obj.N
                
                % evaluate dynamics
                [mu_xkp1,var_xkp1] = obj.f(xk(:,iN),uk(:,iN));
                
                % append dynamics constraints
                ceq_dyn(:,iN+1) = xk(:,iN+1) - mu_xkp1;
                
                % append provided equality constraints(h==0)
                ceq_h(:,iN) = obj.h(xk(:,iN),uk(:,iN));
                
                % provided inequality constraints (g<=0)
                cineq_g(:,iN) = obj.g(xk(:,iN),uk(:,iN));

                t = t + iN * obj.dt;
            end
            ceq   = [ceq_dyn(:); ceq_h(:)];
            cineq = cineq_g(:);
        end
        
    end
end


