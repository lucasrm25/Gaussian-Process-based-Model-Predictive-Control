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
        isOptimized = false % <bool> true if model had its kernel parameters optimized and no new data has been added
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
            assert( obj.n == size(M,1),     'Matrix M has wrong dimension. Expected <%d,%d,%d>',obj.n,obj.n,obj.p);
            assert( obj.p == size(var_f,1), 'Matrix var_f has wrong dimension. Expected <%d>, got <%d>.',obj.p,size(var_f,1));
            assert( obj.p == size(var_n,1), 'Matrix var_n has wrong dimension. Expected <%d>, got <%d>.',obj.p,size(var_n,1));
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
                K=obj.K(obj.X,obj.X);
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
            assert(size(Y,2) == obj.p, ...
                sprintf('Y should have %d columns, but has %d. Dimension does not agree with the specified kernel parameters',obj.p,size(Y,2)));
            assert(size(X,1) == obj.n, ...
                sprintf('X should have %d rows, but has %d. Dimension does not agree with the specified kernel parameters',obj.n,size(X,1)));
            
            
            % dictionary is full
            if obj.N + size(X,2) > obj.Nmax
                % For now, we just keep the most recent data
                % obj.X = [obj.X(:,2:end), X];     % concatenation in the 2st dim.
                % obj.Y = [obj.Y(2:end,:); Y];    % concatenation in the 1st dim.
                
                D = pdist2(obj.X',X','mahalanobis', eye(5) ).^2;
                [~,idx] = max(D);
                
                obj.X = [obj.X(:,1:obj.N ~= idx), X];     % concatenation in the 2st dim.
                obj.Y = [obj.Y(1:obj.N ~= idx,:); Y];    % concatenation in the 1st dim.
                
            % append to dictionary
            else
                obj.X = [obj.X, X];     % concatenation in the 2st dim.
                obj.Y = [obj.Y; Y];     % concatenation in the 1st dim.
            end
        end
        
        
        function [mu_y, var_y] = eval(obj,x)
        %------------------------------------------------------------------
        % Evaluate GP at the points x
        % This is a fast implementation of [Rasmussen, pg19]
        %
        % args:
        %   x: <n,Nx> point coordinates
        % out:
        %   muy:  <p,Nx>      E[Y] = E[gp(x)]
        %   vary: <p,p,Nx>    Var[Y]  = Var[gp(x)]
        %------------------------------------------------------------------
            Nx = size(x,2);  % size of dataset to be evaluated
        
            % if there is no data in the dictionary, return GP prior
            if obj.N == 0 || ~obj.isActive
                mu_y  = repmat(obj.mu(x),[1,obj.p])';
                var_y = zeros(obj.p,obj.p,Nx);
                return;
            end
            
            assert(size(x,1)==obj.n, sprintf('Input vector has %d columns but should have %d !!!',size(x,1),obj.n));
            assert(~isempty(obj.alpha), 'Please call updateModel() at least once before evaluating!!!')
        
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
        
        % function eval_gradx(obj)
        %     KxX = obj.K(x,obj.X);
        %     muy  = obj.mu(x) + KxX * obj.inv_KXX_sn * (obj.Y-obj.mu(obj.X));
        %     vary = obj.K(x,x) - KxX * obj.inv_KXX_sn * KxX';
        % end
        
        
        function optimizeHyperParams(obj)
        %------------------------------------------------------------------
        % Optimize kernel hyper-parameters based on the current dictionary
        %------------------------------------------------------------------
            error('not yet implemented!!!');
            if ~obj.isOptimized
                for ip = 1:obj.p
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % TODO:
                    %       - Implement ML/MAP optimization of hyper parameters
                    %       - See Rasmussen's book Sec. 5.4.1
                    %
                    %       - Each output dimension is a separate GP and must 
                    %       be optimized separately.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                obj.isOptimized = true;
            end
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

