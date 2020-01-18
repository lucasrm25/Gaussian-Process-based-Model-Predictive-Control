
classdef InvertedPendulumGPAnimation < handle
    
    properties       
        % GP
        d_GP
        d_GP_video

        % figure parameters
        x1min
        x1max
        x2min
        x2max
       
        % figure handles
        f
        f_d
        f_sample_points
        f_angle
        
        % mesh parameters
        X1
        X2
        Ytrue
        Ystd
        Ymean
        
        % output simulation
        out
        out_true
        out_nom
        time
        angle
        
        k   % current time step
        
    end
    
    methods
        function obj = InvertedPendulumGPAnimation(d_GP, x1min, x1max, x2min, x2max, out, out_true, out_nom)
         obj.d_GP = d_GP;
         obj.x1min = x1min;
         obj.x1max = x1max;
         obj.x2min = x2min;
         obj.x2max = x2max;
         obj.out = out;
         obj.out_true = out_true;
         obj.out_nom = out_nom;
        end
        
        function initInvertedPendulumGPAnimation(obj)
        % -----------------------------------------------------------------
        %   Generate grid where the mean and variance will be calculated
        % -----------------------------------------------------------------
             % generate double of d_GP
             obj.d_GP_video = GP(obj.d_GP.n, obj.d_GP.p, obj.d_GP.var_f, obj.d_GP.var_n, obj.d_GP.M, obj.d_GP.Nmax);
            %obj.d_GP_video = obj.out_nom.d_GP;

             factor = 0.3;
             rangeX1 = [obj.x1min - factor*range([obj.x1min obj.x1max]), ...
                           obj.x1max + factor*range([obj.x1min obj.x1max])];
             rangeX2 = [obj.x2min - factor*range([obj.x2min obj.x2max]), ...
                           obj.x2max + factor*range([obj.x2min obj.x2max])];

            % generate grid
            [obj.X1,obj.X2] = meshgrid(linspace(rangeX1(1),rangeX1(2),100),...
                               linspace(rangeX2(1),rangeX2(2),100));
            obj.Ytrue = zeros(length(obj.X1),'like',obj.X1);
            obj.Ystd  = zeros(length(obj.X1),'like',obj.X1);
            obj.Ymean = zeros(length(obj.X1),'like',obj.X1);
            
            % figure coordinate system
            obj.f = figure('Color','w','visible','on','units','normalized','outerposition',[0 0 1 1]); grid on; 
            subplot(2,3,[1,2,4,5])
            plot3(NaN,NaN,NaN)
            xlim([min(rangeX1),max(rangeX1)]);
            ylim([min(rangeX2),max(rangeX2)]);
            zlim([-0.01 0.1]);
            hold on
            grid on
            xlabel('$\theta$','Interpreter','latex'); ylabel('$\dot{\theta}$','Interpreter','latex'); zlabel('$\mu(d)$','Interpreter','latex')
            title('Prediction Mean')
            view(190,65)

            
             % -------------------------------------------------------------
            %   setup state predictions
            % -------------------------------------------------------------
            k = 1;
            pi = 1;
            
            for i=1:size(obj.X1,1)
                for j=1:size(obj.X1,2)
                    % evaluate GP model
                    [mu,var] = obj.d_GP_video.eval([obj.X1(i,j);obj.X2(i,j)]);
                    if var < 0
                        error('GP obtained a negative variance... aborting');
                    end
                    obj.Ystd(i,j)  = sqrt(var);
                    obj.Ymean(i,j) = mu(:,pi);    % select desired output dim
                end
            end 
            
            % plot d_GP
            obj.f_d = surf(obj.X1,obj.X2,obj.Ymean, 'FaceAlpha',.8,...
                'EdgeColor', 'none', 'DisplayName', 'Prediction mean',...
                'XDataSource', 'obj.X1',...
                'YDataSource', 'obj.X2',...
                'ZDataSource', 'obj.Ymean');
            hold on
            
            % plot sample points
            obj.f_sample_points = scatter3(NaN,NaN,NaN,...
                     'filled','MarkerFaceColor','red', 'DisplayName', 'Sample points',...
                     'XDataSource', 'obj.d_GP_video.X(1,:)',...
                     'YDataSource', 'obj.d_GP_video.X(2,:)',...
                     'ZDataSource', 'obj.d_GP_video.Y(:,pi)')
            
            % set legend    
            legend([obj.f_d  obj.f_sample_points])
            
            % plot angle
            subplot(2,3,3)
            plot(obj.out_true.out.t(:),obj.out_true.out.xhat(3,:),'DisplayName', 'true model')
            hold on
            plot(obj.out_nom.out.t(:),obj.out_nom.out.xhat(3,:),'DisplayName', 'nominal model')
            hold on
            plot(obj.out.t(1:end-1),obj.out.r(1,:),'DisplayName', 'reference')
            xlim([0,max(obj.out.t)]);
            ylim([min(obj.out.xhat(3,:)),max(obj.out.xhat(3,:))]);
            title('Angle \theta')
            hold on
            grid on
            
            obj.time = repmat(NaN,[1,length(obj.out.t)]);
            obj.angle = repmat(NaN,[1,length(obj.out.t)]);
            obj.f_angle = plot(obj.time, obj.angle,...
                     'DisplayName', '\theta',...
                     'XDataSource', 'obj.time',...
                     'YDataSource', 'obj.angle')
            xlabel('time $t$ [s]','Interpreter','latex');
            ylabel('$\theta$ [rad]','Interpreter','latex');
                     
            legend
  
                     
      
           
            
            % lock up axis limits
            xlim('manual')
            ylim('manual')
            zlim('manual')
            
            drawnow
        end
        
      
        
        
        function status = updateInvertedPendulumGPAnimation(obj,k)
        % -----------------------------------------------------------------
        %   Update GP animation
        % -----------------------------------------------------------------
            pi = 1;
            status = 0;
            if k < 1 
                return;
            end
             
            % add new data points to GP
            obj.d_GP_video.add(obj.d_GP.X(:,k), obj.d_GP.Y(k,:));
            %d_est = [0 0 1 0]' \ (obj.out.xhat(:,k+1) - obj.out.xnom(:,k+1));
            %obj.d_GP_video.add([obj.out.xhat(3,k);obj.out.xhat(4,k)], d_est);
            
            
            for i=1:size(obj.X1,1)
                for j=1:size(obj.X1,2)
                    % evaluate GP model
                    [mu,var] = obj.d_GP_video.eval([obj.X1(i,j);obj.X2(i,j)]);
                    if var < 0
                        error('GP obtained a negative variance... aborting');
                    end
                    obj.Ystd(i,j)  = sqrt(var);
                    obj.Ymean(i,j) = mu(:,pi);    % select desired output dim
                end
            end 
           
           % update
           obj.f_sample_points.XData = obj.d_GP_video.X(1,:);
           obj.f_sample_points.YData = obj.d_GP_video.X(2,:);
           obj.f_sample_points.ZData = obj.d_GP_video.Y(:,pi);
           
           obj.f_d.XData = obj.X1;
           obj.f_d.YData = obj.X2;
           obj.f_d.ZData = obj.Ymean;
           
           obj.f_angle.XData = obj.d_GP_video.X(1,:);
           obj.f_angle.YData = obj.d_GP_video.X(2,:);
           
           obj.time(k) = obj.out.t(k);
           obj.angle(k) = obj.out.xhat(3,k);
           obj.f_angle.XData = obj.time;
           obj.f_angle.YData = obj.angle;
         
           
           
%             obj.f
%             surf(obj.X1,obj.X2,obj.Ymean, 'FaceAlpha',.8, 'EdgeColor', 'none', 'DisplayName', 'Prediction mean');
%             hold on
%             if ~isempty(obj.d_GP_video.X)
%             scatter3(obj.d_GP_video.X(1,:),obj.d_GP_video.X(2,:),obj.d_GP_video.Y(:,pi),'filled','MarkerFaceColor','red', 'DisplayName', 'Sample points')
%             end
%             hold off
            
            % no error when updating graphics
            status = 1;
        end
        
        
        
        function recordvideo(obj, videoName, format, FrameRate)
        % -----------------------------------------------------------------
        %   Record video of the track animation
        % -----------------------------------------------------------------
            % video rec
            videoframes = struct('cdata',[],'colormap',[]);
            obj.initInvertedPendulumGPAnimation();
            xlim manual
            ylim manual
            for k=1:obj.d_GP.N % only holds for N<Nmax
                status = obj.updateInvertedPendulumGPAnimation(k);
                if status == 0
                    break;
                end
                videoframes(k) = getframe(obj.f);
            end
            % -----------------------------------------------------------------
            %   Save video
            % -----------------------------------------------------------------
            writerObj = VideoWriter(videoName,format);
            writerObj.FrameRate = FrameRate;
            open(writerObj);
            % Write out all the frames.
            numberOfFrames = length(videoframes);
            for k=1:numberOfFrames 
               writeVideo(writerObj, videoframes(k));
            end
            close(writerObj);
            disp('Video saved successfully')
        end
    end
end