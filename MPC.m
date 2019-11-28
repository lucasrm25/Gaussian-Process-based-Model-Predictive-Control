classdef MPC
    %------------------------------------------------------------------
    % Nonlinear MPC class
    % 
    %   solve:
    %
    %   min fo(xk,uk,ek)
    %   st  xk+1 = f(xk,uk) + Bd*(d(zk)+wk),    wk~N(0,sigmaw^2)
    %       h(xk,uk,ek) == 0
    %       g(xk,uk,ek) <= 0
    %
    %   for x1,...,xN, u0,...,uN-1, e1,...,eN
    %
    %   where xk: state variables
    %         zk: selected state variables zk=Bd'*xk
    %         ek: extra variables
    %------------------------------------------------------------------
    
    properties
        % Optimizer settings
        maxiter = 200   % max solver iterations
        tol = 1e-8      % optimizer tolerance
        N = 40          % prediction horizon
        
        % Define optimization problem
        fo      % @fun nonlinear cost function
        f       % @fun nonlinear dynamics
        d @GP   % [GP,GP,...] vector of Gaussian process models representing each disturbance dimension
        Bd
        sigmaw
    end
    
    methods
        
        function obj = MPC(f0, f, d, Bd, sigmaw)
        %------------------------------------------------------------------
        % MPC constructor
        % args:
        %   f0: <1> signal/output stddev
        %   f:  <1> evaluation noise stddev
        %   d:  <n,n> length scale covariance matrix
        %   Bd: <1> maximum dictionary size
        %------------------------------------------------------------------
           obj.f0 = f0;
           obj.f  = f;
           obj.d  = d;
           obj.Bd = Bd;
           obj.sigmaw = sigmaw;
        end
        
        
        function u0 = optimize(obj, x0, e0)
            z0 = Bd'*x0;
            u0 = 0;
        end
        
    end
end