%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -
%------------------------------------------------------------------

classdef GP < handle
    %------------------------------------------------------------------
    % Gaussian Process model fitting with SEQ kernel function
    %
    % This is a MIMO fitting function    x \in R^n -> y \in R^p.
    % each output dimension is treated as a different GP
    % 
    % The Kernel parameters, when optimized, are optimized for each output
    % dimension separately.
    %
    %------------------------------------------------------------------
    
    properties
        % kernel parameters (one for each output dimension)
        var_f % <p>     signal/output covariance
        var_n % <p>     evaluation noise covariance
        M     % <n,n,p> length scale covariance
        
        isActive = true
    end
    
    properties(SetAccess=private)
        % dictionary: [X,Y]
        X     % <n,N> input dataset
        Y     % <N,p> output dataset
        
        N     % <1> dictionary size
        Nmax  % <1> max dictionary size
        
        n % <1> input dimension
        p % <1> output dimension
        
        % aux variables
        L       % <N,N,p>: chol(obj.KXX + sigman^2 * I,'lower');
        alpha   % <N,p>: L'\(L\(Y-muX));
        % (DEPRECATED) inv_KXX_sn  % <N,N>
        
        isOutdated  = false % <bool> model is outdated if data has been added withouth updating L and alpha matrices
        % isOptimized = false % <bool> true if model had its kernel parameters optimized and no new data has been added
    end
    
    methods
        function obj = GP(n, p, var_f, var_n, M, maxsize)
        %------------------------------------------------------------------
        % GP constructor
        % args:
        %   n:       <1>     input dimension
        %   p:       <1>     output dimension
        %   var_f:   <p>     signal/output covariance
        %   var_n:   <p>     evaluation noise covariance
        %   M:       <n,n,p> length scale covariance matrix
        %   maxsize: <1>     maximum dictionary size
        %------------------------------------------------------------------
            obj.X       = [];
            obj.Y       = [];
            obj.var_f   = var_f;
            obj.var_n   = var_n;
            obj.M       = M;
            obj.Nmax    = maxsize;
            obj.n       = n;
            obj.p       = p;
            
            % validade model parameters
            assert( obj.n == size(M,1),     'Matrix M has wrong dimension or parameters n/p are wrong. Expected dim(M)=<n,n,p>=<%d,%d,%d>',obj.n,obj.n,obj.p);
            assert( obj.p == size(var_f,1), 'Matrix var_f has wrong dimension or parameter p is wrong. Expected dim(var_f)=<p>=<%d>, got <%d>.',obj.p,size(var_f,1));
            assert( obj.p == size(var_n,1), 'Matrix var_n has wrong dimension or parameter p is wrong. Expected dim(var_n)=<p>=<%d>, got <%d>.',obj.p,size(var_n,1));
        end
        
        
        function bool = isfull(obj)
        %------------------------------------------------------------------
        % is dictionary full?
        %------------------------------------------------------------------
            bool = obj.N >= obj.Nmax;
        end
        
        
        function N = get.N(obj)
        %------------------------------------------------------------------
        % return dictionary size = N
        %------------------------------------------------------------------
            N = size(obj.X,2);
        end
        
        
        function mean = mu(~,x)
        %------------------------------------------------------------------
        % zero mean function: mean[f(x)] = 0
        % args:
        %   x: <n,N>
        % out:
        %   mean: <N,1>
        %------------------------------------------------------------------
            mean = zeros(size(x,2),1);
        end
        
        
        function kernel = K(obj,x1,x2)
        %------------------------------------------------------------------
        % SEQ kernel function: cov[f(x1),f(x2)]
        %     k(x1,x2) = var_f * exp( 0.5 * ||x1-x2||^2_M )
        %
        % args:
        %   x1: <n,N1>
        %   x2: <n,N2>
        % out:
        %   kernel: <N1,N2,p>
        %------------------------------------------------------------------
            for pi = 1:obj.p
                D(:,:,pi) = pdist2(x1',x2','mahalanobis',obj.M(:,:,pi)).^2;
                %D = pdist2(x1',x2','seuclidean',diag((obj.M).^0.5)).^2;
                kernel(:,:,pi) = obj.var_f(pi) * exp( -0.5 * D(:,:,pi) );
            end
        end
        
        
        function updateModel(obj)
        % ----------------------------------------------------------------- 
        % Update precomputed matrices L and alpha, that will be used when
        % evaluating new points. See [Rasmussen, pg19].
        % -----------------------------------------------------------------
            if obj.isOutdated
                % store cholesky L and alpha matrices
                I = eye(obj.N);
                
                % for each output dimension
                obj.alpha = zeros(obj.N,obj.p);
                obj.L = zeros(obj.N,obj.N);
                K = obj.K(obj.X,obj.X);
                for pi=1:obj.p
                    obj.L(:,:,pi) = chol(K(:,:,pi)+ obj.var_n(pi) * I ,'lower');
                    % sanity check: norm( L*L' - (obj.K(obj.X,obj.X) + obj.var_n*I) ) < 1e-12
                    
                    obj.alpha(:,pi) = obj.L(:,:,pi)'\(obj.L(:,:,pi)\(obj.Y(:,pi)-obj.mu(obj.X)));
                end
                
                %-------------------- (DEPRECATED) ------------------------ 
                % % SLOW BUT RETURNS THE FULL COVARIANCE MATRIX INSTEAD OF ONLY THE DIAGONAL (VAR)
                % % precompute inv(K(X,X) + sigman^2*I)
                % I = eye(obj.N);
                % obj.inv_KXX_sn = inv( obj.K(obj.X,obj.X) + obj.var_n * I );
                %-------------------- (DEPRECATED) ------------------------
                
                % set flag
                obj.isOutdated = false;
            end
        end
        
        
        function set.X(obj,X)
            obj.X = X;
            % data has been added. GP is outdated. Please call obj.updateModel
            obj.isOutdated = true;
        end
        
        
        function set.Y(obj,Y)
            obj.Y = Y;
            % data has been added. GP is outdated. Please call obj.updateModel
            obj.isOutdated = true;
        end
        
        
        function add(obj,X,Y)
        %------------------------------------------------------------------
        % Add new data points [X,Y] to the dictionary
        %
        % args:
        %   X: <n,N>
        %   Y: <N,p>
        %------------------------------------------------------------------
            OPTION = 'A'; % {'A','B'}
        
            assert(size(Y,2) == obj.p, ...
                sprintf('Y should have %d columns, but has %d. Dimension does not agree with the specified kernel parameters',obj.p,size(Y,2)));
            assert(size(X,1) == obj.n, ...
                sprintf('X should have %d rows, but has %d. Dimension does not agree with the specified kernel parameters',obj.n,size(X,1)));
            
            Ntoadd = size(X,2);
            Nextra = obj.N + Ntoadd - obj.Nmax;

            % if there is space enough to append the new data points, then
            if Nextra <= 0 
                obj.X = [obj.X, X];
                obj.Y = [obj.Y; Y];
                obj.updateModel();
            
            % data overflow: dictionary will be full. we need to select
            % relevant points
            else 
                
                Nthatfit = obj.Nmax - obj.N;
                
                % make dictionary full
                obj.X = [obj.X, X(:,1:Nthatfit) ];
                obj.Y = [obj.Y; Y(1:Nthatfit,:) ];
                obj.updateModel();
                
                % points left to be added
                X = X(:,Nthatfit+1:end);
                Y = Y(Nthatfit+1:end,:);
                
                % OPTION A) 
                % The closest (euclidian dist.) points will be iteratively removed
                if strcmp(OPTION,'A')
                    for i=1:Nextra
                        D = pdist2(obj.X',X(:,i)','euclidean').^2;
                        [~,idx_rm] = min(D);
                        idx_keep = 1:obj.N ~= idx_rm;

                        obj.X = [obj.X(:,idx_keep), X(:,i)];  % concatenation in the 2st dim.
                        obj.Y = [obj.Y(idx_keep,:); Y(i,:)];    % concatenation in the 1st dim.
                    end
                    
                % OPTION B)
                % the point with lowest variance will be removed        
                else
                    X_all = [obj.X,X];
                    Y_all = [obj.Y;Y];
                    [~, var_y] = obj.eval( X_all, 'activate');
                    [~,idx_keep] = maxk(sum(reshape(var_y, obj.p^2, obj.N+Nextra )),obj.Nmax);

                    obj.X = X_all(:,idx_keep);
                    obj.Y = Y_all(idx_keep,:);
                end
            end
        end
        
        
        function [mu_y, var_y] = eval(obj, x, varargin)
        %------------------------------------------------------------------
        % Evaluate GP at the points x
        % This is a fast implementation of [Rasmussen, pg19]
        %
        % args:
        %   x: <n,Nx> point coordinates
        % varargin: 
        %   'activate': force calculation of mean and variance even if GP is inactive
        % out:
        %   muy:  <p,Nx>      E[Y] = E[gp(x)]
        %   vary: <p,p,Nx>    Var[Y]  = Var[gp(x)]
        %------------------------------------------------------------------
            assert(size(x,1)==obj.n, sprintf('Input vector has %d columns but should have %d !!!',size(x,1),obj.n));
            
            % calculate mean and variance even if GP is inactive
            forceActive = length(varargin)>1 && strcmp(varargin{1},'activate');
            
            Nx = size(x,2);  % size of dataset to be evaluated
        
            % if there is no data in the dictionary or GP is not active, return zeros
            if obj.N == 0 || (~obj.isActive && ~forceActive)
                mu_y  = repmat(obj.mu(x),[1,obj.p])';
                var_y = zeros(obj.p,obj.p,Nx);
                return;
            end
            
            % in case the matrices alpha and L are empty we need to update the model
            if isempty(obj.alpha) || isempty(obj.L)
                obj.updateModel();
            end
            
            % Calculate posterior mean mu_y for each output dimension
            KxX = obj.K(x,obj.X);
            mu_y = zeros(obj.p,Nx);
            for pi=1:obj.p
                 mu_y(pi,:) = obj.mu(x) + KxX(:,:,pi) * obj.alpha(:,pi);
            end
            
            % Calculate posterior covariance var_y
            var_y = zeros(obj.p,obj.p,Nx);
            for pi=1:obj.p
                for i=1:Nx
                    % (less efficient) v = obj.L\obj.K(x(:,i),obj.X)';
                    v = obj.L(:,:,pi)\KxX(i,:,pi)';
                    K = obj.K(x(:,i),x(:,i));
                    var_y(pi,pi,i) = K(:,:,pi) - v'*v;
                end
            end
            
            % --------------------- (DEPRECATED) ------------------------- 
            % % SLOW BUT RETURNS THE FULL COVARIANCE MATRIX INSTEAD OF ONLY THE DIAGONAL (VAR)
            % KxX = obj.K(x,obj.X);
            % muy  = obj.mu(x) + KxX * obj.inv_KXX_sn * (obj.Y-obj.mu(obj.X));
            % vary = obj.K(x,x) - KxX * obj.inv_KXX_sn * KxX';
            % --------------------- (DEPRECATED) ------------------------- 
        end
        
        
        function optimizeHyperParams(obj, method)
        %------------------------------------------------------------------
        % Optimize kernel hyper-parameters based on the current dictionary
        %------------------------------------------------------------------
            
            warning('off', 'MATLAB:nearlySingularMatrix')
            warning('off', 'MATLAB:singularMatrix')
        
            % error('not yet implemented!!!');
            for ip = 1:obj.p
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % TODO:
                %       - Implement ML/MAP optimization of hyper parameters
                %       - See Rasmussen's book Sec. 5.4.1
                %
                %       - Each output dimension is a separate GP and must 
                %       be optimized separately.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % d_GP.optimizeHyperParams
                % obj.optimizeHyperParams_costfun( [obj.var_f; diag(obj.M)])

                nvars = 1 + obj.n;
                IntCon = [];
                fun = @(vars) optimizeHyperParams_costfun(obj,ip,vars);
                A = [];
                b = [];
                Aeq = [];
                beq = [];
                lb = [ 1e-4 ; ones(obj.n,1)*1e-4 ];
                ub = [ 1e+4 ; ones(obj.n,1)*1e+4 ];
                nonlcon = [];

                if strcmp(method, 'ga')
                    options = optimoptions('ga',...
                                           'ConstraintTolerance',1e-6,...
                                           'PlotFcn', @gaplotbestf,...
                                           'Display','iter',...
                                           'UseParallel',false);
                    opt_vars = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options);
                elseif strcmp(method,'fmincon')
                    x0 = [obj.var_f(ip); diag(obj.M(:,:,ip))]; 
                    options = optimoptions('fmincon', ...
                                           'PlotFcn','optimplotfval',...
                                           'Display','iter');
                    [opt_vars,~] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
                else
                    error('Method %s not implemented, please choose an existing method',method);
                end

                obj.var_f(ip) = opt_vars(1);
                obj.M(:,:,ip) = diag( opt_vars(2:end) );


                % update matrices alpha and L
                obj.isOutdated = true;
                obj.updateModel();
            end
            
            warning('on', 'MATLAB:nearlySingularMatrix')
            warning('on', 'MATLAB:singularMatrix')
        end
        
        function cost = optimizeHyperParams_costfun(obj,outdim,vars)
            var_f = vars(1);
            M = diag(vars(2:end));
            Y = obj.Y(:,outdim);
            var_n = obj.var_n(outdim);
            
            K  = var_f * exp( -0.5 * pdist2(obj.X',obj.X','mahalanobis',M).^2 );
            Ky = K + var_n*eye(obj.N);
            
            cost = -0.5* Y' * Ky * Y -0.5* logdet(Ky) - obj.n/2 * log(2*pi);
        end
        
        function cost = optimizeHyperParams_gradfun(obj,outdim,vars)
            
            var_f = vars(1);
            M = diag(vars(2:end));
                
            K = var_f * exp( -0.5 * pdist2(obj.X',obj.X','mahalanobis',M).^2 );
            
            alpha = K \ obj.Y(:,outdim);
            
            dK_var_f = K*2/sqrt(var_f);
            
            dK_l = zeros(obj.N,obj.N);
            for i=1:obj.N
                for j=1:obj.N
                    ksi = obj.X(:,i) - obj.X(:,j);
                    % dK_l(i,j) = sum( K(i,j)*0.5*inv(M)*ksi*ksi'*inv(M) * log(diag(M)) );
                    dK_l(i,j) = sum( K(i,j)*0.5*M\ksi*ksi'/M * log(diag(M)) );
                end
            end            
            % cost = 0.5 * trace( (alpha*alpha' - inv(K)) * ( dK_var_f + dK_l ) );
            cost = 0.5 * trace( alpha*alpha'*(dK_var_f+dK_l) - K\(dK_var_f+dK_l) );
        end
        
        
        function plot2d(obj, truthfun, varargin)
        %------------------------------------------------------------------
        % Make analysis of the GP quality (only for the first output dimension.
        % This function can only be called when the GP input is 2D
        %
        % args:
        %   truthfun: anonymous function @(x) which returns the true function
        %   varargin{1} = rangeX1: 
        %   varargin{2} = rangeX2:  <1,2> range of X1 and X2 where the data 
        %                           will be evaluated and ploted
        %------------------------------------------------------------------
            % output dimension to be analyzed
            pi = 1;
        
            assert(obj.N>0, 'Dataset is empty. Aborting...')
            % we can not plot more than in 3D
            assert(obj.n==2, 'This function can only be used when dim(X)=2. Aborting...');
            
            % Generate grid where the mean and variance will be calculated
            if numel(varargin) ~= 2
                factor = 0.3;
                rangeX1 = [ min(obj.X(1,:)) - factor*range(obj.X(1,:)), ...
                            max(obj.X(1,:)) + factor*range(obj.X(1,:))  ];
                rangeX2 = [ min(obj.X(2,:)) - factor*range(obj.X(2,:)), ...
                            max(obj.X(2,:)) + factor*range(obj.X(2,:))  ];
            else
                rangeX1 = varargin{1};
                rangeX2 = varargin{2};
            end

            % generate grid
            [X1,X2] = meshgrid(linspace(rangeX1(1),rangeX1(2),100),...
                               linspace(rangeX2(1),rangeX2(2),100));
            Ytrue = zeros('like',X1);
            Ystd  = zeros('like',X1);
            Ymean = zeros('like',X1);
            for i=1:size(X1,1)
                for j=1:size(X1,2)
                    % evaluate true function
                    mutrue = truthfun([X1(i,j);X2(i,j)]);
                    Ytrue(i,j) = mutrue(pi); % select desired output dim
                    % evaluate GP model
                    [mu,var] = obj.eval([X1(i,j);X2(i,j)]);
                    if var < 0
                        error('GP obtained a negative variance... aborting');
                    end
                    Ystd(i,j)  = sqrt(var);
                    Ymean(i,j) = mu(:,pi);    % select desired output dim
                end
            end 
            
            % plot data points, and +-2*stddev surfaces 
            figure('Color','w', 'Position', [-1827 27 550 420])
            hold on; grid on;
            % surf(X1,X2,Y, 'FaceAlpha',0.3)
            surf(X1,X2,Ymean+2*Ystd ,Ystd, 'FaceAlpha',0.3)
            surf(X1,X2,Ymean-2*Ystd,Ystd, 'FaceAlpha',0.3)
            scatter3(obj.X(1,:),obj.X(2,:),obj.Y(:,pi),'filled','MarkerFaceColor','red')
            title('mean\pm2*stddev Prediction Curves')
            shading interp;
            colormap(gcf,jet);
            view(30,30)
            
            % Comparison between true and prediction mean
            figure('Color','w', 'Position',[-1269 32 1148 423])
            subplot(1,2,1); hold on; grid on;
            surf(X1,X2,Ytrue, 'FaceAlpha',.8, 'EdgeColor', 'none', 'DisplayName', 'True function');
            % surf(X1,X2,Ymean, 'FaceAlpha',.5, 'FaceColor','g', 'EdgeColor', 'none', 'DisplayName', 'Prediction mean');
            scatter3(obj.X(1,:),obj.X(2,:),obj.Y(:,pi),'filled','MarkerFaceColor','red', 'DisplayName', 'Sample points')
            zlim([ min(obj.Y(:,pi))-range(obj.Y(:,pi)),max(obj.Y(:,pi))+range(obj.Y(:,pi)) ]);
            legend;
            xlabel('X1'); ylabel('X2');
            title('True Function')
            view(24,12)
            subplot(1,2,2); hold on; grid on;
            % surf(X1,X2,Y, 'FaceAlpha',.5, 'FaceColor','b', 'EdgeColor', 'none', 'DisplayName', 'True function');
            surf(X1,X2,Ymean, 'FaceAlpha',.8, 'EdgeColor', 'none', 'DisplayName', 'Prediction mean');
            scatter3(obj.X(1,:),obj.X(2,:),obj.Y(:,pi),'filled','MarkerFaceColor','red', 'DisplayName', 'Sample points')
            zlim([ min(obj.Y(:,pi))-range(obj.Y(:,pi)),max(obj.Y(:,pi))+range(obj.Y(:,pi)) ]);
            legend;
            xlabel('X1'); ylabel('X2');
            title('Prediction Mean')
            view(24,12)
            
            % plot bias and variance
            figure('Color','w', 'Position',[-1260 547 894 264])
            subplot(1,2,1); hold on; grid on;
            contourf(X1,X2, abs(Ymean-Ytrue), 50,'LineColor','none')
            title('Absolute Prediction Bias')
            colorbar;
            scatter(obj.X(1,:),obj.X(2,:),'filled','MarkerFaceColor','red')
            subplot(1,2,2); hold on; grid on;
            contourf(X1,X2, Ystd.^2, 50 ,'LineColor','none')
            title('Prediction Variance')
            colorbar;
            scatter(obj.X(1,:),obj.X(2,:),'filled','MarkerFaceColor','red')
            colormap(gcf,parula);
        end
        
        
        function plot1d(obj, truthfun, varargin)
        %------------------------------------------------------------------
        % Make analysis of the GP quality (only for the first output dimension.
        % This function can only be called when the GP input is 1D
        %
        % args:
        %   truthfun: anonymous function @(x) which returns the true function
        %   varargin{1} = rangeX1: 
        %   varargin{2} = rangeX2:  <1,2> range of X1 and X2 where the data 
        %                           will be evaluated and ploted
        %------------------------------------------------------------------
        
            assert(obj.N>0, 'Dataset is empty. Aborting...')
            % we can not plot more than in 3D
            assert(obj.n==1, 'This function can only be used when dim(X)=1. Aborting...');
            
            error('Not implemented error')
        end
    end
end

