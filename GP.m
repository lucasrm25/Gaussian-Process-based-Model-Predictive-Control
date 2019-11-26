classdef GP < handle
    
    properties
        dict
        sigmaf
        sigman
        L
        maxsize
        inv_KXX_sn
    end
    
    methods
        function obj = GP(sigmaf, sigman, lambda, maxsize)
        %------------------------------------------------------------------
        % GP constructor
        % args:
        %   sigmaf: <1> signal/output stddev
        %   sigman: <1> evaluation noise stddev
        %   lambda: <n,n> length scale covariance matrix
        %   maxsize: maximum dictionary size
        %------------------------------------------------------------------
            obj.dict.X  = [];
            obj.dict.Y  = [];
            obj.sigmaf  = sigmaf;
            obj.sigman  = sigman;
            obj.L       = lambda;
            obj.maxsize = maxsize;
        end
        
        
        function mean = mu(~,x)
        %------------------------------------------------------------------
        % zero mean function: mean[f(x)] = 0
        % args:
        %   x: <D,N>
        %------------------------------------------------------------------
            mean = zeros(size(x,2),1);
        end
        
        
        function kernel = K(obj,x1,x2)
        %------------------------------------------------------------------
        % SEQ kernel function: cov[f(x1),f(x2)]
        % args:
        %   x1: <D,N1>
        %   x2: <D,N2>
        % out:
        %   kernel: <N2,N1>
        %------------------------------------------------------------------
            nx1 = size(x1,2);
            nx2 = size(x2,2);
            kernel = zeros(nx1,nx2);
            for i=1:nx1
                for j=1:nx2
                    % kernel(i,j) = obj.sigmaf^2 * exp( - norm(x1(:,i)-x2(:,j))^2 / (2*obj.L^2) );
                    kernel(i,j) = obj.sigmaf^2 * exp( -0.5 * (x1(:,i)-x2(:,j))'/obj.L*(x1(:,i)-x2(:,j)) );
                end
            end
        end
        
        
        function add(obj,X,Y)
        %------------------------------------------------------------------
        % args:
        %   X: <D,N>
        %   Y: <1,N>
        %------------------------------------------------------------------
            if size(obj.dict.X,2) + size(X,2) > obj.maxsize
                % decide which points to select
            else
                obj.dict.X = [obj.dict.X, X];
                obj.dict.Y = [obj.dict.Y; Y];
                % precompute inv(K(X,X) + sigman^2*I)
                I = eye(size(obj.dict.X,2));
                obj.inv_KXX_sn = inv( obj.K(obj.dict.X,obj.dict.X) + obj.sigman^2 * I );
            end
            
        end
        
        
        function [muy, covary] = eval(obj,x)
        %------------------------------------------------------------------
        % evaluate GP at the points x
        % args:
        %   x: <D,N> point coordinates
        %------------------------------------------------------------------
            KxX = obj.K(x,obj.dict.X);
            KXx = KxX';
            muy  = obj.mu(x) + KxX * obj.inv_KXX_sn * (obj.dict.Y-obj.mu(obj.dict.X));
            covary = obj.K(x,x) - KxX * obj.inv_KXX_sn * KXx;
        end
        
        
        function plot2d(obj, truthfun, varargin)
        %------------------------------------------------------------------
        % plot mean and covariance of GP with 2D states
        % args:
        %   truthfun: anonymous function @(x) which returns the true function
        %   varargin{1} = rangeX1, 
        %   varargin{2} = rangeX2:  <1,2> range of X1 and X2 where the data 
        %                           will be evaluated
        %------------------------------------------------------------------   
            if size(obj.dict.X,1) ~= 2
                error('This function can only be used when dim(X)=2');
            end
            
            
            % Generate grid where the mean and variance will be calculated
            if numel(varargin) ~= 2
                factor = 0.3;
                rangeX1 = [ min(obj.dict.X(1,:)) - factor*range(obj.dict.X(1,:)), ...
                            max(obj.dict.X(1,:)) + factor*range(obj.dict.X(1,:))  ];
                rangeX2 = [ min(obj.dict.X(2,:)) - factor*range(obj.dict.X(2,:)), ...
                            max(obj.dict.X(2,:)) + factor*range(obj.dict.X(2,:))  ];
            else
                rangeX1 = varargin{1};
                rangeX2 = varargin{2};
            end

            % generate grid
            [X1,X2] = meshgrid(linspace(rangeX1(1),rangeX1(2),100),...
                               linspace(rangeX2(1),rangeX2(2),100));
            for i=1:size(X1,1)
                for j=1:size(X1,2)
                    % evaluate true function
                    Y(i,j) = truthfun([X1(i,j);X2(i,j)]);
                    % evaluate GP model
                    [mu,var] = obj.eval([X1(i,j);X2(i,j)]);
                    if var < 0
                        error('GP obtained a negative variance... aborting');
                    end
                    Ystd(i,j)  = sqrt(var);
                    Ymean(i,j) = mu;
                end
            end 
            
            
            % plot data points, and +-2*stddev surfaces 
            figure('Color','w', 'Position', [123   124   550   420])
            hold on; grid on;
            % surf(X1,X2,Y, 'FaceAlpha',0.3)
            surf(X1,X2,Ymean+2*Ystd ,Ystd, 'FaceAlpha',0.3)
            surf(X1,X2,Ymean-2*Ystd,Ystd, 'FaceAlpha',0.3)
            scatter3(obj.dict.X(1,:),obj.dict.X(2,:),obj.dict.Y,'filled','MarkerFaceColor','red')
            title('mean\pm2*stddev Prediction Curves')
            shading interp;
            colormap(gcf,jet);
            view(30,30)
            
            % Comparison between true and prediction mean
            figure('Color','w', 'Position',[286 146 1138 423])
            subplot(1,2,1); hold on; grid on;
            surf(X1,X2,Y, 'FaceAlpha',.8, 'EdgeColor', 'none', 'DisplayName', 'True function');
            % surf(X1,X2,Ymean, 'FaceAlpha',.5, 'FaceColor','g', 'EdgeColor', 'none', 'DisplayName', 'Prediction mean');
            scatter3(obj.dict.X(1,:),obj.dict.X(2,:),obj.dict.Y,'filled','MarkerFaceColor','red', 'DisplayName', 'Sample points')
            zlim([ min(obj.dict.Y)-range(obj.dict.Y),max(obj.dict.Y)+range(obj.dict.Y) ]);
            legend;
            xlabel('X1'); ylabel('X2');
            title('True Function')
            view(24,12)
            subplot(1,2,2); hold on; grid on;
            % surf(X1,X2,Y, 'FaceAlpha',.5, 'FaceColor','b', 'EdgeColor', 'none', 'DisplayName', 'True function');
            surf(X1,X2,Ymean, 'FaceAlpha',.8, 'EdgeColor', 'none', 'DisplayName', 'Prediction mean');
            scatter3(obj.dict.X(1,:),obj.dict.X(2,:),obj.dict.Y,'filled','MarkerFaceColor','red', 'DisplayName', 'Sample points')
            zlim([ min(obj.dict.Y)-range(obj.dict.Y),max(obj.dict.Y)+range(obj.dict.Y) ]);
            legend;
            xlabel('X1'); ylabel('X2');
            title('Prediction Mean')
            view(24,12)
            
            % plot bias and variance
            figure('Color','w', 'Position',[708   166   894   264])
            subplot(1,2,1); hold on; grid on;
            contourf(X1,X2, abs(Ymean-Y), 50,'LineColor','none')
            title('Absolute Prediction Bias')
            colorbar;
            scatter(obj.dict.X(1,:),obj.dict.X(2,:),'filled','MarkerFaceColor','red')
            subplot(1,2,2); hold on; grid on;
            contourf(X1,X2, Ystd.^2, 50 ,'LineColor','none')
            title('Prediction Variance')
            colorbar;
            scatter(obj.dict.X(1,:),obj.dict.X(2,:),'filled','MarkerFaceColor','red')
            colormap(gcf,parula);
        end
    end
end

