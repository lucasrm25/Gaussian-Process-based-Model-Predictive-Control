%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -
%------------------------------------------------------------------

classdef NMPC
    %------------------------------------------------------------------
    % Nonlinear MPC class
    % 
    %   solve:
    %
    %   MIN     { SUM_{k=0:N-1} fo(tk,xk,uk,ek,r(tk)) } + fend(tN,xN,eN,r(tN))
    %
    %   s.t.    xk+1 = f(xk,uk) + Bd*(d(zk)+wk),    wk~N(0,sigmaw^2)
    %           h(tk,xk,uk,ek) == 0
    %           g(tk,xk,uk,ek) <= 0
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
        tol     = 1e-8      % optimizer tolerance
        N       = 30        % prediction horizon
        
        % Define optimization problem
        fo      % @fun nonlinear cost function
        fend    % @fend nonlinear cost function for the final state
        f       % @fun nonlinear dynamics
        d @GP   % [GP,GP,...] vector of Gaussian process models representing each disturbance dimension
        Bd
        sigmaw  % discrete noise covariance of w
        h       % equality constraint function 
        g       % inequality constraint function
        
        % Optimization dimension
        n   % state space dimension
        m   % input diension
        ne  % extra variable dimension
        
        % time step
        dt  % time step size
        
        isGPactive = false  % if inactive, NMPC will consider d(zk)=0
    end
    
    properties(Access=private)
        % save last optimal results computed, in order to use as initial guess
        vars_opt_old = []
    end
    
    methods
        
        function obj = NMPC(fo, fend, f, d, Bd, N, sigmaw, h, g, n, m, ne, dt)
        %------------------------------------------------------------------
        % MPC constructor
        % args:
        %   f0: <1> signal/output stddev
        %   f:  <1> evaluation noise stddev
        %   d:  <n,n> length scale covariance matrix
        %   Bd: <1> maximum dictionary size
        %------------------------------------------------------------------
           obj.fo = fo;
           obj.f  = f;
           obj.fend = fend;
           obj.d  = d;
           obj.Bd = Bd;
           obj.N = N;
           obj.sigmaw = sigmaw;
           obj.h = h;
           obj.g = g;
           obj.n = n;
           obj.m = m;
           obj.ne = ne;
           obj.dt = dt;
        end
        
        
        function u0_opt = optimize(obj, x0, e0, t0, r)
            z0 = obj.Bd'*x0;
            
            if isempty(obj.vars_opt_old)
                % initialize optimization variables initial guesses
                xguess = x0*ones(obj.n * obj.N+1,1);     % [ x0,...,xN-1,xN ]
                uguess = zeros(obj.m * obj.N,1);         % [ u0,...,uN-1    ]
                eguess = e0*ones(obj.ne * obj.N+1,1);    % [ e0,...,eN-1,eN ]
                % vector of initial guesses
                varsguess = [xguess; uguess; eguess];
            else
                % use last optimization results as initial guess in case
                % they have already been calculated once
                varsguess = obj.vars_opt_old;
            end
                
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
                                    'Display','none',...
                                   'Algorithm','sqp',...
                                   'ConstraintTolerance',obj.tol,...
                                   'MaxIterations',obj.maxiter);
            
            % solve optimization problem                   
            vars_opt = fmincon(costfun,varsguess,A,b,Aeq,beq,lb,ub,nonlcon,options);
            
            % store current optimization results to use as initial guess
            % for future optimizations
            obj.vars_opt_old = vars_opt;
            
            % split variables since vars_opt = [x_opt; u_opt; e_opt]
            [x_opt, u_opt, e_opt] = splitvariables(obj, vars_opt);
            u0_opt = u_opt(1:obj.m);
            
            % plot results every time step
            if false
                figure; hold on; grid on;
                plot(x_opt, 'DisplayName','x');
                plot(u_opt, 'DisplayName','u');
            end
        end
        
    
        function [xvec, uvec, evec] = splitvariables(obj, vars)
            % split variables
            xvec = vars(1:(obj.N+1)*obj.n);
            uvec = vars( (1:obj.N*obj.m)  + length(xvec) );
            evec = vars( (1:(obj.N+1)*obj.ne) + length(xvec) + length(uvec) );
        end
        

        function cost = costfun(obj, vars, t0, r)
            % to make the notation shorter
            n = obj.n;
            m = obj.m;
            ne = obj.ne;
            N = obj.N;

            % split variables
            [xvec, uvec, evec] = obj.splitvariables(vars);

            cost = 0;
            t = t0;
            for iN=1:N+1      % i=0:N
                
                % which are the indices of the current time k?
                idx_xk =  n*(iN-1)+1 : n*iN;
                idx_uk =  m*(iN-1)+1 : m*iN;
                idx_ek = ne*(iN-1)+1 : ne*iN;

                % calculate cost for t = tk
                if iN <= N
                    % cost for current state
                    cost = cost + obj.fo( t, xvec(idx_xk), uvec(idx_uk), evec(idx_ek), r );
                else
                    % final state cost
                    cost = cost + obj.fend( t, xvec(idx_xk), evec(idx_ek), r );
                end
                
                % update current time
                t = t + iN * obj.dt;
            end
            
        end
        

        function [cineq,ceq] = nonlcon(obj,vars, t0, x0)

            % to make the notation shorter
            n = obj.n;
            m = obj.m;
            ne = obj.ne;
            N = obj.N;

            % split variables
            [xvec, uvec, evec] = obj.splitvariables(vars);
            
            % initial state constraint: x0 - x(0) = 0
            ceq = [ x0 - xvec(1) ];
            cineq = [];
            
            t = t0;
            for iN=1:N

                % which are the indices of the current time k?
                idx_xk   =  n*(iN-1)+1 : n*iN;
                idx_uk   =  m*(iN-1)+1 : m*iN;
                idx_ek   = ne*(iN-1)+1 : ne*iN;

                
                if obj.isGPactive
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % TODO:
                    %       calculate constraint when GP is active
                    %       see Formula Student AMZ paper
                    %       mean and covariance needs to be propagated with
                    %       EKF
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                else % add constraints without GP disturbance estimation
                    
                    % dynamics (xk+1-f(xk,uk)=0) and provided equality (h=0) constraints
                    ceq = [ceq ;
                           xvec(idx_xk+n) - obj.f(t,xvec(idx_xk),uvec(idx_uk)) ;
                           obj.h(t,xvec(idx_xk),uvec(idx_uk),evec(idx_ek))     ];
                    
                    % provided inequality constraints (g<=0)   
                    cineq = [cineq;
                             obj.g(t,xvec(idx_xk),uvec(idx_uk),evec(idx_ek)) ];
                end
                t = t + iN * obj.dt;
            end

        end
        
    end
end


