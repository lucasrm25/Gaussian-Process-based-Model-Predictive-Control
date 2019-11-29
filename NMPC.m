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
        tol     = 1e-6      % optimizer tolerance
        N       = 30        % prediction horizon
        
        % Define optimization problem
        fo      % @fun nonlinear cost function
        fend    % @fend nonlinear cost function for the final state
        f       % @fun nonlinear dynamics
        d @GP   % Gaussian Process model representing disturbance
        Bd
        sigmaw  % discrete noise covariance of w
        h       % equality constraint function 
        g       % inequality constraint function
        
        % Optimization dimension
        n   % dim(x) state space dimension
        m   % dim(u) input diension
        ne  % dim(e) extra variable dimension
        nz  % dim(z) selected variable dimension - z = Bd*x
        dt  % time step size
    end
    
    properties(Access=private)
        % save last optimal results computed, in order to use as initial guess
        vars_opt_old = []
        
        isGPactive = false  % if inactive, NMPC will consider d(zk)=0
    end
    
    methods
        
        function obj = NMPC(fo, fend, f, d, Bd, N, sigmaw, h, g, n, m, ne, dt)
        %------------------------------------------------------------------
        % MPC constructor
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
           obj.nz = size(Bd,1);
        end
        
        function activateGP(obj)
            if obj.d.N == 0
                error('GP dataset is empty. Please add data points before activating!');
            end
            obj.isGPactive = true;
        end
        
        function deactivateGP(obj)
            obj.isGPactive = false;
        end
        
        
        function numvars = optSize(obj)
        %------------------------------------------------------------------
        % How many variables we need to optimize?
        %
        %   vars_opt = [x0,...,xN, u0,...,uN-1, e0,...,eN]
        %------------------------------------------------------------------
            numvars = (obj.N+1)*(obj.n) + (obj.N)*(obj.m) + (obj.N+1)*(obj.ne);
        end
        
        
        function u0_opt = optimize(obj, x0, e0, t0, r)
        %------------------------------------------------------------------
        % Calculate first uptimal control input
        %------------------------------------------------------------------
            
        
            %-------- Set initial guess for optimization variables  -------
            % initialize optimization variables initial guesses
            if isempty(obj.vars_opt_old)
                % if this is the first optimization
                xguess = x0*ones(obj.n * (obj.N+1),1);     % [ x0,...,xN-1,xN ]
                uguess = zeros(obj.m * obj.N,1);           % [ u0,...,uN-1    ]
                eguess = e0*ones(obj.ne * (obj.N+1),1);    % [ e0,...,eN-1,eN ]
                % vector of initial guesses
                varsguess = [xguess; uguess; eguess];
            else
                varsguess = obj.vars_opt_old;
                varsguess(1:obj.n) = x0;
                
                % use last optimization results as initial guess in case
                % they have already been calculated once
                % % [xvec, uvec, evec] = obj.splitvariables(obj.vars_opt_old);
                % % 
                % % varsguess = [ x0; xvec(2*obj.n+1:end); xvec(end-obj.n+1:end);
                % %               uvec(obj.m+1:end); uvec(end-obj.m+1:end)];
                % % if ~isempty(evec) && numel(evec) > 1
                % %     varsguess = [varguess; evec(2:end); evec(end)];
                % % else
                % %     varsguess = [varsguess; evec];
                % % end
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
            
            % if ~isempty(obj.vars_opt_old)
            %     [xerr, uerr, eerr] = obj.splitvariables(vars_opt-obj.vars_opt_old);
            %     fprintf('\nOld guess error. x:%.3f - u:%.3f - e: %.3f\n', ...
            %          norm(xerr), norm(uerr), norm(eerr) );
            %      [xerr, uerr, eerr] = obj.splitvariables(vars_opt-varsguess);
            %     fprintf('\nNew guess error. x:%.3f - u:%.3f - e: %.3f\n', ...
            %          norm(xerr), norm(uerr), norm(eerr) );
            % end
            
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
        %------------------------------------------------------------------
        % Evaluate cost function for the whole horizon, given variables
        %------------------------------------------------------------------

            % split variables
            [xvec, uvec, evec] = obj.splitvariables(vars);

            cost = 0;
            t = t0;
            for iN=1:obj.N+1      % i=0:N
                
                % which are the indices of the current vector at time k?
                idx_xk =  obj.n*(iN-1)+1 : obj.n*iN;
                idx_uk =  obj.m*(iN-1)+1 : obj.m*iN;
                idx_ek = obj.ne*(iN-1)+1 : obj.ne*iN;

                % calculate cost for t = tk
                if iN <= obj.N
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
        %------------------------------------------------------------------
        % Evaluate nonlinear equality and inequality constraints
        % args
        %   cineq = g(x,u)
        %   ceq   = h(x,u)
        %------------------------------------------------------------------ 

            % split variables
            [xvec, uvec, evec] = obj.splitvariables(vars);
            
            % initialize state constraint: x0 - x(0) = 0
            ceq = x0 - xvec(1);
            % initialize inequality constraint
            cineq = [];
            
            if obj.isGPactive
                % evaluate GP for all points in the horizon
                zvec = obj.Bd*xvec;
                [mu_dz,var_dz] = obj.d.eval([zvec(1:end-obj.nz),uvec]');
            end
            
            t = t0;
            for iN=1:obj.N

                % which are the indices of the current time k?
                idx_xk   = obj.n *(iN-1)+1 : obj.n *iN;
                idx_zk   = obj.nz*(iN-1)+1 : obj.nz*iN;
                idx_uk   = obj.m *(iN-1)+1 : obj.m *iN;
                idx_ek   = obj.ne*(iN-1)+1 : obj.ne*iN;
                
                % get current timestep variables
                xk   = xvec(idx_xk);
                xkp1 = xvec(idx_xk + obj.n);
                uk   = uvec(idx_uk);
                ek   = evec(idx_ek);
                zk   = zvec(idx_zk); %obj.Bd * xk;
                
                if obj.isGPactive
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % TODO:
                    %       calculate constraint when GP is active
                    %       see Formula Student AMZ paper
                    %       mean and covariance needs to be propagated with
                    %       EKF
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % dynamics constraints
                    % (DEPRECATED) [mu_dzk,var_dzk] = obj.d.eval([zk;uk]);
                    mu_dzk = mu_dz(idx_zk);
                    ceq = [ceq ;
                           xkp1 - obj.f(t,xk,uk) - obj.Bd * mu_dzk];
                       
                    % cineq = [cineq;
                    %         ];
                    
                else % add constraints without GP disturbance estimation
                    
                    % dynamics constraints
                    ceq = [ceq ;
                           xkp1 - obj.f(t,xk,uk) ];
                end
                
                % provided equality constraints (h==0)
                if ~isempty(obj.h)
                    ceq = [ceq ; obj.h(t,xk,uk,ek) ];
                end
                
                % provided inequality constraints (g<=0)
                if ~isempty(obj.g)
                    cineq = [cineq; obj.g(t,xk,uk,ek) ];
                end
                
                t = t + iN * obj.dt;
            end

        end
        
    end
end


