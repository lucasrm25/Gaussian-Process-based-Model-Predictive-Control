%------------------------------------------------------------------
% Programed by: 
%   - Lucas Rath (lucasrm25@gmail.com)
%   - 
%   -
%------------------------------------------------------------------

classdef GP < handle
    %------------------------------------------------------------------
    % Gaussian Process model fitting with SEQ kernel function
    %------------------------------------------------------------------
    
    
    properties
        % dictionary: [X,Y]
        X        % <D,N>
        Y        % <1,N>
        maxsize  % <1>
        
        % kernel parameters
        sigmaf  % <1>
        sigman  % <1>
        L       % <D,D>
    end
    
    properties(Access=private)
        % aux variables
        inv_KXX_sn  % <N,N>
    end
    
    methods
        function obj = GP(sigmaf, sigman, lambda, maxsize)
        %------------------------------------------------------------------
        % GP constructor
        % args:
        %   sigmaf: <1> signal/output stddev
        %   sigman: <1> evaluation noise stddev
        %   lambda: <n,n> length scale covariance matrix
        %   maxsize: <1> maximum dictionary size
        %------------------------------------------------------------------
            obj.X  = [];
            obj.Y  = [];
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
        %   kernel: <N1,N2>
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
        
        function set.X(obj, X)
        %------------------------------------------------------------------
        % precompute inv(K(X,X) + sigman^2*I), every time X changes
        %------------------------------------------------------------------
            % set variable
            obj.X = X;
            % precompute inv(K(X,X) + sigman^2*I)
            I = eye(size(obj.X,2));
            obj.inv_KXX_sn = inv( obj.K(obj.X,obj.X) + obj.sigman^2 * I );
        end
        
        
        function add(obj,X,Y)
        %------------------------------------------------------------------
        % args:
        %   X: <D,N>
        %   Y: <1,N>
        %------------------------------------------------------------------
            if size(obj.X,2) + size(X,2) > obj.maxsize
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % TODO:
                %       - decide how to select the induction points
                %       - READ the paper from AMZ. They give hints there
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                obj.X = [obj.X, X];
                obj.Y = [obj.Y; Y];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TODO:
            %       - call optimizeHyperParams ???
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        
        function [muy, covary] = eval(obj,x)
        %------------------------------------------------------------------
        % evaluate GP at the points x
        % args:
        %   x: <D,N> point coordinates
        %------------------------------------------------------------------
            KxX = obj.K(x,obj.X);
            muy  = obj.mu(x) + KxX * obj.inv_KXX_sn * (obj.Y-obj.mu(obj.X));
            covary = obj.K(x,x) - KxX * obj.inv_KXX_sn * KxX';
        end
        
        
        function optimizeHyperParams(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TODO:
            %       - Implement ML/MAP optimization of hyper parameters
            %       - See Rasmussen's book Sec. 5.4.1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            return;
        end
        
        
        function plot2d(obj, truthfun, varargin)
        %------------------------------------------------------------------
        % plot mean and covariance of GP with 2D states
        % args:
        %   truthfun: anonymous function @(x) which returns the true function
        %   varargin{1} = rangeX1: 
        %   varargin{2} = rangeX2:  <1,2> range of X1 and X2 where the data 
        %                           will be evaluated and ploted
        %------------------------------------------------------------------   
            if size(obj.X,1) ~= 2
                error('This function can only be used when dim(X)=2');
            end
            
            
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
                    Ytrue(i,j) = truthfun([X1(i,j);X2(i,j)]);
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
            scatter3(obj.X(1,:),obj.X(2,:),obj.Y,'filled','MarkerFaceColor','red')
            title('mean\pm2*stddev Prediction Curves')
            shading interp;
            colormap(gcf,jet);
            view(30,30)
            
            % Comparison between true and prediction mean
            figure('Color','w', 'Position',[286 146 1138 423])
            subplot(1,2,1); hold on; grid on;
            surf(X1,X2,Ytrue, 'FaceAlpha',.8, 'EdgeColor', 'none', 'DisplayName', 'True function');
            % surf(X1,X2,Ymean, 'FaceAlpha',.5, 'FaceColor','g', 'EdgeColor', 'none', 'DisplayName', 'Prediction mean');
            scatter3(obj.X(1,:),obj.X(2,:),obj.Y,'filled','MarkerFaceColor','red', 'DisplayName', 'Sample points')
            zlim([ min(obj.Y)-range(obj.Y),max(obj.Y)+range(obj.Y) ]);
            legend;
            xlabel('X1'); ylabel('X2');
            title('True Function')
            view(24,12)
            subplot(1,2,2); hold on; grid on;
            % surf(X1,X2,Y, 'FaceAlpha',.5, 'FaceColor','b', 'EdgeColor', 'none', 'DisplayName', 'True function');
            surf(X1,X2,Ymean, 'FaceAlpha',.8, 'EdgeColor', 'none', 'DisplayName', 'Prediction mean');
            scatter3(obj.X(1,:),obj.X(2,:),obj.Y,'filled','MarkerFaceColor','red', 'DisplayName', 'Sample points')
            zlim([ min(obj.Y)-range(obj.Y),max(obj.Y)+range(obj.Y) ]);
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

