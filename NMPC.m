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
    %   for x0, u0,...,uN-1, e1,...,eN
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
        h       % nonlinear equality constraint function 
        g       % nonlinear inequality constraint function        
        u_lb    % inputs lower bound constraint  u_lb <= u
        u_ub    % inputs upper bound constraint          u <= u_ub
        
        fo      % @fun nonlinear cost function
        fend    % @fend nonlinear cost function for the final state
        
        % Optimization dimension
        n   % dim(x) state space dimension
        m   % dim(u) input dimension
        ne  % dim(e) additional optimization variables dimension
        dt  % time step size
        nh  % number of additional eq. constraints for every time step
        ng  % number of additional ineq. constraints for every time step
    end
    
    properties(Access=private)
        % save last optimal results computed, in order to use as initial guess
        vars_opt_old = []
    end
    
    
    methods
        
        function obj = NMPC (f, h, g, u_lb, u_ub, n, m, ne, fo, fend, N, dt)
        %------------------------------------------------------------------
        % MPC constructor
        %
        %   f: motion model that evaluates  [E[xk+1],Var[xk+1]] = f(xk)
        %------------------------------------------------------------------
           % constraints
           obj.f  = f;
           obj.h = h;
           obj.g = g;
           obj.u_lb = u_lb;
           obj.u_ub = u_ub;
           % variable dimensions
           obj.n = n;
           obj.m = m;
           obj.ne = ne;
           % get size of additional constraints
           obj.nh = length(h(zeros(n,1),zeros(m,1),zeros(ne,1)));
           obj.ng = length(g(zeros(n,1),zeros(m,1),zeros(ne,1)));
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
            numvars = obj.n + obj.N*obj.m + obj.N*obj.ne;
        end
        
        
        function [x0_opt, u_opt, e_opt] = optimize(obj, x0, t0, r)
        %------------------------------------------------------------------
        % Calculate first uptimal control input
        %------------------------------------------------------------------
            
            %-------- Set initial guess for optimization variables  -------
            % initialize optimization variables initial guesses
            if isempty(obj.vars_opt_old)  % if this is the first optimization
                uguess = zeros(obj.m*obj.N, 1);      % [ u0,...,uN-1]
                eguess = zeros(obj.ne*obj.N, 1);     % [ e1,...,eN]
            else
                [~,uguess,eguess] = obj.splitvariables(obj.vars_opt_old);
                uguess = uguess(:); % <m,N>  to <m*N,1>
                eguess = eguess(:); % <ne,N> to <ne*N,1>
            end
            varsguess = [x0; uguess; eguess];
            %--------------------------------------------------------------
            
            assert( numel(varsguess) == obj.optSize(), ...
                'There is something wrong with the code. Number of optimization variables does not match!' );
                
            % define cost and constr. functions, as a function only of the
            % optimazation variables. This is a prerequisite for the function fmincon
            costfun = @(vars) obj.costfun(vars, t0, r);
            nonlcon = @(vars) obj.nonlcon(vars, t0, x0);
            
            % define lower and upper bound constraints
            lb = [-Inf(obj.n,1);
                   repmat(obj.u_lb,obj.N,1);    % repeat lower bound for all u0,...,uN-1
                   -Inf(obj.ne*obj.N,1)];             
               
            ub = [Inf(obj.n,1);
                  repmat(obj.u_ub,obj.N,1);     % repeat upper bound for all u0,...,uN-1
                  Inf(obj.ne*obj.N,1)];
            
            % define optimizer settings
            options = optimoptions('fmincon',...
                                   'Display','iter',...
                                   'Algorithm','interior-point',... % 'spq'
                                   'ConstraintTolerance',obj.tol,...
                                   'MaxIterations',obj.maxiter);
            
            % solve optimization problem                   
            vars_opt = fmincon(costfun,varsguess,[],[],[],[],lb,ub,nonlcon,options);

            % store current optimization results to use as initial guess
            % for future optimizations
            obj.vars_opt_old = vars_opt;
            
            % split variables since vars_opt = [x_opt; u_opt; e_opt]
            [x0_opt, u_opt, e_opt] = splitvariables(obj, vars_opt);
            % u0_opt = u_opt(:,1);
        end
        
    
        function [x0, uvec, evec] = splitvariables(obj, vars)
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
            evec = vars( (1:obj.N*obj.ne) + length(x0)+ length(uvec) );
            % reshape the column vector <m*N,1> to <m,N>
            uvec = reshape(uvec, obj.m, obj.N);
            % reshape the column vector <ne*N,1> to <ne,N>
            evec = reshape(evec, obj.ne, obj.N);
        end
        
        
        function [mu_xk,var_xk] = predictStateSequence(obj, mu_x0, var_x0, uk)
        %------------------------------------------------------------------
        % Propagate mean and covariance of state sequence, given control
        % input sequence.
        % out:
        %   mu_xk:  <n,N+1>
        %   var_xk: <n,n,N+1>
        %------------------------------------------------------------------
            mu_xk  = zeros(obj.n,obj.N+1);
            var_xk = zeros(obj.n,obj.n,obj.N+1);
            
            mu_xk(:,1) = mu_x0;
            var_xk(:,:,1) = var_x0;
            
            for iN=1:obj.N      % [x1,...,xN]
                try
                    [mu_xk(:,iN+1),var_xk(:,:,iN+1)] = obj.f(mu_xk(:,iN),var_xk(:,:,iN),uk(:,iN));
                catch e
                    error('System dynamics evaluated to error!!!')
                end
                if sum(isnan(mu_xk),'all') || sum(isinf(mu_xk),'all')
                    error('System dynamics evaluated to NaN or Inf')
                end
            end
        end

        
        function cost = costfun(obj, vars, t0, r)
        %------------------------------------------------------------------
        % Evaluate cost function for the whole horizon, given variables
        %------------------------------------------------------------------
            % split variables
            [mu_x0, uvec, evec] = obj.splitvariables(vars);
            var_x0 = zeros(obj.n);
            
            % calculate state sequence for given control input sequence and x0
            [mu_xvec,var_xvec] = obj.predictStateSequence(mu_x0, var_x0, uvec);
                        
            cost = 0;
            t = t0;
            for iN=1:obj.N      % i=0:N-1
                % add cost: fo=@(t,mu_x,var_x,u,e,r)
                try
                    cost = cost + obj.fo(t, mu_xvec(:,iN), var_xvec(:,:,iN), uvec(:,iN), evec(:,iN), r);
                catch e
                    error('Cost function evaluated to error!!!')
                end
                if sum(isnan(cost),'all') || sum(isinf(cost),'all')
                    error('Cost function evaluated to NaN or Inf')
                end
                % update current time
                t = t + iN * obj.dt;
            end
            % final cost: fend=@(t,mu_x,var_x,e,r)
            cost = cost + obj.fend(t, mu_xvec(:,end), var_xvec(:,:,end), evec(:,iN), r);
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
            [mu_x0, uvec, evec] = obj.splitvariables(vars);
            var_x0 = zeros(obj.n);
            
            % calculate state sequence for given control input sequence and x0
            [mu_xvec,var_xvec] = obj.predictStateSequence(mu_x0, var_x0, uvec);
            
            % set initial state constraint: x0 - x(0) = 0
            ceq_dyn(:,1) = x0 - mu_xvec(:,1);
            
            t = t0;
            for iN=1:obj.N
                try
                    % append provided equality constraints(h==0)
                    ceq_h(:,iN) = obj.h(mu_xvec(:,iN),uvec(:,iN), evec(:,iN));
                    % provided inequality constraints (g<=0)
                    cineq_g(:,iN) = obj.g(mu_xvec(:,iN),uvec(:,iN),evec(:,iN));
                catch e
                    error('Constraints h(h) or g(x) evaluated to error!!!')
                end

                t = t + iN * obj.dt;
            end
            ceq   = [ceq_dyn(:); ceq_h(:)];
            cineq = cineq_g(:);
        end
        
    end
end


