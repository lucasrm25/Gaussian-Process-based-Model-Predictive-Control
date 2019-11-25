clear all; close all; clc;



% Hyperparameters
sigmaf  = 10;        % output variance (std)
lambda  = 10;        % length scale
sigman  = 0.05;     % STD of measurement noise
maxsize = 100;      % maximum number of points in the dictionary

gp = GP(sigmaf, sigman, lambda, maxsize);



truthfun = @(x) x'*[1;2] - 20;
xdata = rand(2,30)*30;
ydata = truthfun(xdata) + sigman*randn(size(xdata,2),1);

gp.add(xdata,ydata)


x = [25:35; 25:35];
[muy, covary] = gp.eval(x);



%% plots

gp.plot2d(-5:40,-5:40)

% [X1,X2] = meshgrid(-5:40,-5:40);
% for i=1:size(X1,1)
%     for j=1:size(X1,2)
%         Y(i,j) = truthfun([X1(i,j);X2(i,j)]);
%         [mu,var] = gp.eval([X1(i,j);X2(i,j)]);
%         Yl(i,j) = mu - 2*sqrt(var);
%         Yu(i,j) = mu + 2*sqrt(var);
%         Ystd(i,j) = sqrt(var);
%     end
% end
% 
% figure('Color','w')
% hold on; grid on;
% % surf(X1,X2,Y, 'FaceAlpha',0.3)
% surf(X1,X2,Yl,Ystd, 'FaceAlpha',0.3)
% surf(X1,X2,Yu,Ystd, 'FaceAlpha',0.3)
% scatter3(xdata(1,:),xdata(2,:),ydata,'filled')
% 
% shading interp;
% colormap jet
% view(30,30)