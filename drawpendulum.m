function drawpendulum(t,x,t2,x2, m1, m2,g,l)
%% drawpendulum
% adapted from Patrick Suhms code 
% (source: https://github.com/PatrickSuhm/LqrControlTutorial)
% ------------------------------------------------------------------------
% input:
% t    ... time (,1)
% x    ... state ([s, ds, phi, dphi])
% m1   ... mass of the cart
% m2   ... mass of pole
% ------------------------------------------------------------------------

s = x(1,:);                             % position   
phi = x(3,:);                           % angle


s2 = x2(1,:);                             % position   
phi2 = x2(3,:);   

% dimensions of cart and mass
W = 1*sqrt(m1/5);  % cart width
H = .5*sqrt(m1/5); % cart height
mr = .3*sqrt(m2);  % mass radius

% position mass
px = s - l*sin(phi);
py = H/2 + l*cos(phi);

px2 = s2 - l*sin(phi2);
py2 = H/2 + l*cos(phi2);

% create new figure and 
figure
plot([-25 25],[0 0],'w','LineWidth',2)
hold on

% plot the cart
h1=rectangle('Position',[s(1)-W/2,0,W,H],'Curvature',.1,'FaceColor',[1 0.1 0.1],'EdgeColor',[1 1 1]);
h12=rectangle('Position',[s2(1)-W/2,0,W,H],'Curvature',.1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);    
    % plot the pole
h2=plot([s(1) px(1)],[0 py(1)],'w','LineWidth',2);
h22=plot([s2(1) px(1)],[0 py2(1)],'w','LineWidth',2);

h3=rectangle('Position',[px(1)-mr/2,py(1)-mr/2,mr,mr],'Curvature',1,'FaceColor',[.3 0.3 1],'EdgeColor',[1 1 1]);
h32=rectangle('Position',[px2(1)-mr/2,py2(1)-mr/2,mr,mr],'Curvature',1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1]);

h4=text(3, 4.5, ['time: ', num2str(1)]);
% h5=text(2, 4, ['angle (nom): ', num2str(1)]);
% h6=text(2, 3.5, ['angle (GP): ', num2str(1)]);

set(h4,'color','w', 'fontsize', 14);
% set(h5,'color','w', 'fontsize', 14);
% set(h6,'color','w', 'fontsize', 14);
xlim([-5 5]);
ylim([-2 5]);
set(gca,'Color','k','XColor','w','YColor','w')
set(gcf,'Color','k')

pause(0.2) 
tic

videoframes = struct('cdata',[],'colormap',[]);
% animation in a for loop 
for k=1:length(t)

  % update pole and cart position
  set(h1, 'position',[s(k)-W/2,0,W,H]);
  set(h2, 'XData',[s(k) px(k)], 'YData', [H/2 py(k)]);
  set(h3, 'position', [px(k)-mr/2,py(k)-mr/2,mr,mr]);
  set(h4, 'string',['time: ', num2str(t(k))]);
%   set(h5, 'string',['angle (nom): ', num2str(rad2deg(phi(k)))]);
%   set(h6, 'string',['angle (GP): ', num2str(rad2deg(phi2(k)))]);
  
  % set second cart
  set(h12, 'position',[s2(k)-W/2,0,W,H]);
  set(h22, 'XData',[s2(k) px2(k)], 'YData', [H/2 py2(k)]);
  set(h32, 'position', [px2(k)-mr/2,py2(k)-mr/2,mr,mr]);
  videoframes(k) = getframe(gcf);
  drawnow();
  pause(0.2)
  
  % strop the animation when q is pressed
%   if k==length(t)
%     stop
%   end
  
  % meassure the time and create a fixed time loop
%   t2=toc;
%   while t2 < t(k)
%     t2 = toc;
%   end
%   t3(k) = t2;
  
end
FrameRate = 5;
videoName = fullfile('simresults',sprintf('InvertedPendulumSim-%s',date));
videoFormat = 'Motion JPEG AVI';

 writerObj = VideoWriter(videoName,videoFormat);
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
