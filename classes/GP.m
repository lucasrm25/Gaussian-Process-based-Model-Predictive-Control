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
        isActive = true
    end
    
    properties(SetAccess=private)
        % kernel parameters (one for each output dimension)
        M     % <n,n,p> length scale covariance
        var_f % <p>     signal/output covariance
        var_n % <p>     evaluation noise covariance
        
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
        
        % isOutdated  = false % <bool> model is outdated if data has been added withouth updating L and alpha matrices
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
            obj.n       = n;
            obj.p       = p;
            obj.X       = [];
            obj.Y       = [];
            obj.Nmax    = maxsize;    
            obj.setHyperParameters( M, var_f, var_n )
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
        
        function setHyperParameters(obj, M, var_f, var_n)
        %------------------------------------------------------------------
        % set new values to the hyperparameter
        %------------------------------------------------------------------    
            assert( obj.n == size(M,1),     'Matrix M has wrong dimension or parameters n/p are wrong. Expected dim(M)=<n,n,p>=<%d,%d,%d>',obj.n,obj.n,obj.p);
            assert( obj.p == size(var_f,1), 'Matrix var_f has wrong dimension or parameter p is wrong. Expected dim(var_f)=<p>=<%d>, got <%d>.',obj.p,size(var_f,1));
            assert( obj.p == size(var_n,1), 'Matrix var_n has wrong dimension or parameter p is wrong. Expected dim(var_n)=<p>=<%d>, got <%d>.',obj.p,size(var_n,1));
            obj.M     = M;
            obj.var_f = var_f;
            obj.var_n = var_n;
            obj.updateModel();
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
        
        
        function kernel = K(obj, x1, x2)
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
            kernel = zeros(size(x1,2),size(x2,2),obj.p);
            for pi = 1:obj.p
                D = pdist2(x1',x2','mahalanobis',obj.M(:,:,pi)).^2;
                %D = pdist2(x1',x2','seuclidean',diag((obj.M).^0.5)).^2;
                kernel(:,:,pi) = obj.var_f(pi) * exp( -0.5 * D );
            end
        end
        
        
        function updateModel(obj)
        % ----------------------------------------------------------------- 
        % Update precomputed matrices L and alpha, that will be used when
        % evaluating new points. See [Rasmussen, pg19].
        % -----------------------------------------------------------------
            % nothing to update... return
            if obj.N == 0
                return;
            end
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
        end
        
        
        function add(obj, X, Y)
        %------------------------------------------------------------------
        % Add new data points [X,Y] to the dictionary
        %
        % args:
        %   X: <n,N>
        %   Y: <N,p>
        %------------------------------------------------------------------
            OPTION = 'B'; % {'A','B'}
        
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

                        obj.X = [obj.X(:,idx_keep), X(:,i)];
                        obj.Y = [obj.Y(idx_keep,:); Y(i,:)];
                    end
                    
                % OPTION B)
                % the point with lowest variance will be removed        
                elseif strcmp(OPTION,'B')
                    X_all = [obj.X,X];
                    Y_all = [obj.Y;Y];
                    
                    [~, var_y] = obj.eval( X_all, true);
                    [~,idx_keep] = maxk(sum(reshape(var_y, obj.p^2, obj.N+Nextra )),obj.Nmax);

                    obj.X = X_all(:,idx_keep);
                    obj.Y = Y_all(idx_keep,:);
                else
                    error('Option not implemented');
                end
            end
            % update pre-computed matrices
            obj.updateModel();
        end
        
        
        function [mu_y, var_y] = eval(obj, x, varargin)
        %------------------------------------------------------------------
        % Evaluate GP at the points x
        % This is a fast implementation of [Rasmussen, pg19]
        %
        % args:
        %   x: <n,Nx> point coordinates
        % varargin: 
        %   true: force calculation of mean and variance even if GP is inactive
        % out:
        %   muy:  <p,Nx>      E[Y] = E[gp(x)]
        %   vary: <p,p,Nx>    Var[Y]  = Var[gp(x)]
        %------------------------------------------------------------------
            assert(size(x,1)==obj.n, sprintf('Input vector has %d columns but should have %d !!!',size(x,1),obj.n));
            
            % calculate mean and variance even if GP is inactive
            forceActive = length(varargin)>=1 && varargin{1}==true;
            
            Nx = size(x,2);  % size of dataset to be evaluated
        
            % if there is no data in the dictionary or GP is not active
            % then return prior (for now returning zero variance)
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
        %   Method: maximum log likelihood (See Rasmussen's book Sec. 5.4.1)
        %------------------------------------------------------------------
            
            warning('off', 'MATLAB:nearlySingularMatrix')
            warning('off', 'MATLAB:singularMatrix')
            
            obj.updateModel();
        
            % error('not yet implemented!!!');
            for ip = 1:obj.p
                
                % set optimization problem
                nvars = obj.n + 1 + 1;  % M, var_f, var_n
                %fun = @(vars) loglikelihood(obj,ip,vars);
                
                fun = @(vars) loglikelihood(obj,ip,...                  % output dim
                                            diag(10.^vars(1:end-2)),... % M
                                            10.^vars(end-1),...         % var_f
                                            10.^vars(end));             % var_n
                
                ub = [ 1e+5 * ones(obj.n,1);
                       1e+5;
                       1e+5 ];
                lb = [ 1e-8 * ones(obj.n,1); 
                       0*1e-8; 
                       1e-10 ];
                x0 = [ diag(obj.M(:,:,ip)); obj.var_f(ip); obj.var_n(ip); ];
                
                % convert to log10 space
                ub = log10(ub);
                lb = log10(lb);
                x0 = log10(x0);

                % use genetic algorithm or interior-point-method
                if strcmp(method, 'ga')
                    options = optimoptions('ga',...
                                           'ConstraintTolerance',1e-6,...
                                           'PlotFcn', @gaplotbestf,...
                                           'Display','iter',...
                                           'UseParallel',false);
                    opt_vars = ga(fun,nvars,[],[],[],[],lb,ub,[],[],options);
                elseif strcmp(method,'fmincon')
                    options = optimoptions('fmincon', ...
                                           'PlotFcn','optimplotfval',...
                                           'Display','iter');
                    [opt_vars,~] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
                else
                    error('Method %s not implemented, please choose an existing method',method);
                end

                % retrieve optimal results to absolute scale
                obj.M(:,:,ip) = diag( 10.^opt_vars(1:end-2) );
                obj.var_f(ip) = 10.^opt_vars(end-1);
                obj.var_n(ip) = 10.^opt_vars(end);
            end
            
            % update matrices alpha and L
            obj.updateModel();
            
            warning('on', 'MATLAB:nearlySingularMatrix')
            warning('on', 'MATLAB:singularMatrix')
        end
        
        function logL = loglikelihood(obj, outdim, M, var_f, var_n)
        %------------------------------------------------------------------
        % calculate the negative log likelihood: -log(p(Y|X,theta)),
        %      where theta are the hyperparameters and (X,Y) the training data
        %------------------------------------------------------------------         
            Y = obj.Y(:,outdim);
            K  = var_f * exp( -0.5 * pdist2(obj.X',obj.X','mahalanobis',M).^2 );
            Ky = K + var_n*eye(obj.N);
            % calculate log(p(Y|X,theta))
            logL = -(-0.5*Y'/Ky*Y -0.5*logdet(Ky) -obj.n/2*log(2*pi));
        end
        
        function plot2d(obj, truefun, varargin)
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
            
            assert(obj.N>0, 'Dataset is empty. Aborting...')
            assert(obj.n==2, 'This function can only be used when dim(X)=2. Aborting...');
            
            %--------------------------------------------------------------
            % parse inputs
            %--------------------------------------------------------------
            p = inputParser;
            
            addParameter(p,'factor',0.3); 
            addParameter(p,'outdim',1);
            addParameter(p,'npoints',50);
            addParameter(p,'sigmalevel',2);
            parse(p,varargin{:});
            
            factor = p.Results.factor;
            outdim = p.Results.outdim;
            npoints = p.Results.npoints;
            sigmalevel = p.Results.sigmalevel;
            
            addParameter(p,'rangeX1', minmax(obj.X(1,:)) + [-1 1]*factor*range(obj.X(1,:)) );
            addParameter(p,'rangeX2', minmax(obj.X(2,:)) + [-1 1]*factor*range(obj.X(2,:)) );
            parse(p,varargin{:});
            
            rangeX1 = p.Results.rangeX1;
            rangeX2 = p.Results.rangeX2;
            
            %--------------------------------------------------------------
            % Evaluate Ytrue, Ymean and Ystd
            %--------------------------------------------------------------
            
            % generate grid
            [X1,X2] = meshgrid(linspace(rangeX1(1),rangeX1(2),npoints),...
                               linspace(rangeX2(1),rangeX2(2),npoints));
            Ytrue = zeros('like',X1);
            Ystd  = zeros('like',X1);
            Ymean = zeros('like',X1);
            for i=1:size(X1,1)
                for j=1:size(X1,2)
                    % evaluate true function
                    mutrue = truefun([X1(i,j);X2(i,j)]);
                    Ytrue(i,j) = mutrue(outdim); % select desired output dim
                    % evaluate GP model
                    [mu,var] = obj.eval([X1(i,j);X2(i,j)],true);
                    if var < 0
                        error('GP obtained a negative variance... aborting');
                    end
                    Ystd(i,j)  = sqrt(var);
                    Ymean(i,j) = mu(:,outdim);    % select desired output dim
                end
            end 
            
            %--------------------------------------------------------------
            % Generate plots
            %--------------------------------------------------------------
            
            % plot data points, and +-2*stddev surfaces 
            figure('Color','w', 'Position', [-1827 27 550 420])
            %figure('Color','white','Position',[513  440  560  420]);
            hold on; grid on;
            s1 = surf(X1, X2, Ymean,'edgecolor',0.8*[1 1 1], 'EdgeAlpha', 0.3 ,'FaceColor', [153, 51, 255]/255);
            s2 = surf(X1, X2, Ymean+sigmalevel*Ystd, Ystd, 'FaceAlpha',0.2,'EdgeAlpha',0.2, 'EdgeColor',0.4*[1 1 1]); %, 'FaceColor',0*[1 1 1])
            s3 = surf(X1, X2, Ymean-sigmalevel*Ystd, Ystd,  'FaceAlpha',0.2,'EdgeAlpha',0.2, 'EdgeColor',0.4*[1 1 1]); %, 'FaceColor',0*[1 1 1]) 
            p1 = scatter3(obj.X(1,:),obj.X(2,:),obj.Y(:,outdim),'filled','MarkerFaceColor','red');
            title('mean\pm2*stddev Prediction Curves')
            xlabel('X1'); ylabel('X2');
            view(70,10)
            colormap(gcf,jet);
            
            
            % Comparison between true and prediction mean
            figure('Color','w', 'Position',[-1269 32 1148 423])
            subplot(1,2,1); hold on; grid on;
            surf(X1,X2,Ytrue, 'FaceAlpha',.8, 'EdgeColor', 'none', 'DisplayName', 'True function');
            % surf(X1,X2,Ymean, 'FaceAlpha',.5, 'FaceColor','g', 'EdgeColor', 'none', 'DisplayName', 'Prediction mean');
            scatter3(obj.X(1,:),obj.X(2,:),obj.Y(:,outdim),'filled','MarkerFaceColor','red', 'DisplayName', 'Sample points')
            zlim( minmax(obj.Y(:,outdim)') +[-1 1]*range(obj.Y(:,outdim)) );
            legend;
            title('True Function')
            xlabel('X1'); ylabel('X2');
            view(-60,17)
            subplot(1,2,2); hold on; grid on;
            % surf(X1,X2,Y, 'FaceAlpha',.5, 'FaceColor','b', 'EdgeColor', 'none', 'DisplayName', 'True function');
            surf(X1,X2,Ymean, 'FaceAlpha',.8, 'EdgeColor', 'none', 'DisplayName', 'Prediction mean');
            scatter3(obj.X(1,:),obj.X(2,:),obj.Y(:,outdim),'filled','MarkerFaceColor','red', 'DisplayName', 'Sample points')
            zlim( minmax(obj.Y(:,outdim)') +[-1 1]*range(obj.Y(:,outdim)) );
            legend;
            title('Prediction Mean')
            xlabel('X1'); ylabel('X2');
            view(-60,17)
            
            % plot bias and variance
            figure('Color','w', 'Position',[-1260 547 894 264])
            subplot(1,2,1); hold on; grid on;
            contourf(X1,X2, abs(Ymean-Ytrue), 50,'LineColor','none')
            title('Absolute Prediction Bias')
            xlabel('X1'); ylabel('X2');
            colorbar;
            scatter(obj.X(1,:),obj.X(2,:),'filled','MarkerFaceColor','red')
            subplot(1,2,2); hold on; grid on;
            contourf(X1,X2, Ystd.^2, 50 ,'LineColor','none')
            title('Prediction Variance')
            xlabel('X1'); ylabel('X2');
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
        %   see inputParser parameters
        %------------------------------------------------------------------
        
            assert(obj.N>0, 'Dataset is empty. Aborting...')
            % we can not plot more than in 3D
            assert(obj.n==1, 'This function can only be used when dim(X)=1. Aborting...');
            
            %--------------------------------------------------------------
            % parse inputs
            %--------------------------------------------------------------
            p = inputParser;
            
            addParameter(p,'factor',0.3); 
            addParameter(p,'outdim',1);
            addParameter(p,'npoints',300);
            addParameter(p,'sigmalevel',2);
            parse(p,varargin{:});
            
            factor = p.Results.factor;
            outdim = p.Results.outdim;
            npoints = p.Results.npoints;
            sigmalevel = p.Results.sigmalevel;
            
            addParameter(p,'rangeX', minmax(obj.X) + [-1 1]*factor*range(obj.X) );
            parse(p,varargin{:});
            
            rangeX = p.Results.rangeX;
            
            
            %--------------------------------------------------------------
            % Evaluate Ytrue, Ymean and Ystd
            %--------------------------------------------------------------
            
            % generate grid
            X = linspace(rangeX(1),rangeX(2),npoints);
            % evaluate and calculate prediction mean+-2*std
            [mu,var] = obj.eval(X,true);
            Ytrue = truthfun(X);
            Ymean = mu';
            Ystd  = sqrt(squeeze(var));
            
            % prior
            %Ymean = 0*mu';
            %Ystd  = sqrt(diag(obj.K(X,X)));
            
            %--------------------------------------------------------------
            % Generate plots
            %--------------------------------------------------------------
            
            figure('Color','w'); hold on; grid on;
            p0 = plot(X,Ytrue,      'LineWidth',2);
            p1 = plot(X,Ymean,      'LineWidth',0.5,'Color', [77, 0, 153]/255);
            p2 = plot(X,Ymean+sigmalevel*Ystd, 'LineWidth',0.5,'Color', [77, 0, 153]/255);
            p3 = plot(X,Ymean-sigmalevel*Ystd, 'LineWidth',0.5,'Color', [77, 0, 153]/255);
            p4 = patch([X fliplr(X)], [Ymean'+sigmalevel*Ystd' fliplr(Ymean'-sigmalevel*Ystd')], [153, 51, 255]/255, ...
                        'FaceAlpha',0.2, 'EdgeColor','none');
            p5 = scatter( obj.X, obj.Y, 'MarkerFaceColor','r','MarkerEdgeColor','r');
            % title('mean \pm 2*std curves');
            xlabel('X'); ylabel('Y');
        end
    end
end
