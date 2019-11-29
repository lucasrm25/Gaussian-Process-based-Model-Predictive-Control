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
    %------------------------------------------------------------------
    
    properties
        % kernel parameters (one for each output dimension)
        sigmaf  % <p> signal/output stddev
        sigman  % <p> evaluation noise stddev
        l       % <n,n,p> length scale covariance matrix
    end
    
    properties(SetAccess=private)
        % dictionary: [X,Y]
        X        % <n,N> input dataset
        Y        % <N,p> output dataset
        
        N     % <1> dictionary size
        Nmax  % <1> max dictionary size
        
        n % <1> input dimension
        p % <1> output dimension
        
        % aux variables
        L       % <N,N>: chol(obj.KXX + sigman^2 * I,'lower');
        alpha   % <N,1>: L'\(L\(Y-muX));
        % (DEPRECATED) inv_KXX_sn  % <N,N>
        
        outdated % <bool> tells if data has been added withouth updating L and alpha matrices
    end
    
    methods
        function obj = GP(sigmaf, sigman, lambda, maxsize)
        %------------------------------------------------------------------
        % GP constructor
        % args:
        %   sigmaf: <p> signal/output stddev
        %   sigman: <p> evaluation noise stddev
        %   lambda: <n,n,p> length scale covariance matrix
        %   maxsize: <1> maximum dictionary size
        %------------------------------------------------------------------
            obj.X       = [];
            obj.Y       = [];
            obj.sigmaf  = sigmaf;
            obj.sigman  = sigman;
            obj.l       = lambda;
            obj.Nmax    = maxsize;
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
        
        function n = get.n(obj)
        %------------------------------------------------------------------
        % return input dimension = n
        %------------------------------------------------------------------
            n = size(obj.X,1);
        end
        
        function p = get.p(obj)
        %------------------------------------------------------------------
        % return output dimension
        %------------------------------------------------------------------
            p = size(obj.Y,2);
        end
        
        function mean = mu(~,x)
        %------------------------------------------------------------------
        % zero mean function: mean[f(x)] = 0
        % args:
        %   x: <n,N>
        %------------------------------------------------------------------
            mean = zeros(size(x,2),1);
        end
        
        
        function kernel = K(obj,x1,x2)
        %------------------------------------------------------------------
        % SEQ kernel function: cov[f(x1),f(x2)]
        % args:
        %   x1: <n,N1>
        %   x2: <n,N2>
        % out:
        %   kernel: <N1,N2>
        %------------------------------------------------------------------
            
            % ---------------- (DEPRECATED - TOO SLOW) --------------------
            % nx1 = size(x1,2);
            % nx2 = size(x2,2);
            % kernel = zeros(nx1,nx2);
            % invl = inv(obj.l);
            % for i=1:nx1
            %     for j=1:nx2
            %         dx = x1(:,i) - x2(:,j);
            %         kernel(i,j) = obj.sigmaf^2 * exp( -0.5 * dx'*invl*dx );
            %     end
            % end
            % ---------------- (DEPRECATED - TOO SLOW) --------------------
            
            % ---------------- (DEPRECATED - TOO SLOW) --------------------
            %D = pdist2(x1',x2','mahalanobis',obj.l).^2;
            %kernel = obj.sigmaf^2 * exp( -0.5 * D );
            % ---------------- (DEPRECATED - TOO SLOW) --------------------
            
            D = pdist2(x1',x2','seuclidean',diag((obj.l).^0.5)).^2;
            kernel = obj.sigmaf^2 * exp( -0.5 * D );
        end
        
        function updateModel(obj)
        % ----------------------------------------------------------------- 
        % Update precomputed matrices L and alpha
        % -----------------------------------------------------------------
            if obj.outdated
                obj.outdated = false;
                % store cholesky L and alpha matrices
                I = eye(obj.N);
                obj.L = chol( obj.K(obj.X,obj.X) + obj.sigman^2 * I ,'lower');
                % sanity check: norm( L*L' - (obj.K(obj.X,obj.X) + obj.sigman^2*I) ) < 1e-12
                obj.alpha = obj.L'\(obj.L\(obj.Y-obj.mu(obj.X)));
                
                %-------------------- (DEPRECATED) ------------------------ 
                % % SLOW BUT RETURNS THE FULL COVARIANCE MATRIX INSTEAD OF ONLY THE DIAGONAL (VAR)
                % % precompute inv(K(X,X) + sigman^2*I)
                % I = eye(obj.N);
                % obj.inv_KXX_sn = inv( obj.K(obj.X,obj.X) + obj.sigman^2 * I );
                %-------------------- (DEPRECATED) ------------------------

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % TODO:
                %       - call optimizeHyperParams every time data is added ???
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
        
        function L = get.L(obj)
            if obj.outdated
                obj.updateModel();
            end
            L = obj.L;
        end
        
        function alpha = get.alpha(obj)
            if obj.outdated
                obj.updateModel();
            end
            alpha = obj.alpha;
        end
        
        function set.X(obj,X)
            obj.X = X;
            % data has been added. GP is outdated. Please call obj.updateModel
            obj.outdated = true;
        end
        
        function set.Y(obj,Y)
            obj.Y = Y;
            % data has been added. GP is outdated. Please call obj.updateModel
            obj.outdated = true;
        end
        
        function add(obj,X,Y)
        %------------------------------------------------------------------
        % - add new data points [X,Y] to the dictionary
        % - precompute cholesky L and alpha matrices that will be used when
        %   evaluating new points. See [Rasmussen, pg19].
        %
        % - (DEPRECATED) precompute inv(K(X,X) + sigman^2*I), every time X changes
        %
        % args:
        %   X: <n,N>
        %   Y: <N,p>
        %------------------------------------------------------------------
            if obj.N + size(X,2) > obj.Nmax
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % TODO:
                %       - decide how to select the induction points
                %       - READ the paper from AMZ. They give hints there
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                obj.X = [obj.X, X];     % concatenation in the 2st dim.
                obj.Y = [obj.Y; Y];     % concatenation in the 1st dim.
            end
        end
        
        
        function [muy, vary] = eval(obj,x)
        %------------------------------------------------------------------
        % Evaluate GP at the points x
        % This is a fast implementation of [Rasmussen, pg19]
        %
        % args:
        %   x: <n,N> point coordinates
        % out:
        %   muy:  <N,p>    E[Y]
        %   vary: <N,1>  Var[Y] is the same for all output dimensions
        %   (DEPRECATED) covary: <N,N>
        %------------------------------------------------------------------
            if obj.N == 0
                error('GP dataset is empty. Please add data points before evaluating!');
            end
        
            muy = obj.mu(x) + obj.K(x,obj.X) * obj.alpha;
            
            Nx = size(x,2);  % size of dataset to be evaluated
            vary = zeros(Nx,1);
            for i=1:Nx
                v = obj.L\obj.K(x(:,i),obj.X)';
                vary(i) = obj.K(x(:,i),x(:,i)) - v'*v;
            end
            
            % --------------------- (DEPRECATED) ------------------------- 
            % % SLOW BUT RETURNS THE FULL COVARIANCE MATRIX INSTEAD OF ONLY THE DIAGONAL (VAR)
            % KxX = obj.K(x,obj.X);
            % muy  = obj.mu(x) + KxX * obj.inv_KXX_sn * (obj.Y-obj.mu(obj.X));
            % covary = obj.K(x,x) - KxX * obj.inv_KXX_sn * KxX';
            % --------------------- (DEPRECATED) ------------------------- 
        end
        
        
        function optimizeHyperParams(obj)
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
            return;
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
            figure('Color','w', 'Position', [123   124   550   420])
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
            figure('Color','w', 'Position',[286 146 1138 423])
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
            figure('Color','w', 'Position',[708   166   894   264])
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
    end
end

