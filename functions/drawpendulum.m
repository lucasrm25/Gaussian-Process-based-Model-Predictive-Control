function drawpendulum(t,x, m1, m2,g,l)
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


% dimensions of cart and mass
W = 1*sqrt(m1/5);  % cart width
H = .5*sqrt(m1/5); % cart height
mr = .3*sqrt(m2);  % mass radius

% position mass
px = s - l*sin(phi);
py = H/2 + l*cos(phi);

% create new figure and 
figure
plot([-25 25],[0 0],'w','LineWidth',2)
hold on

% plot the cart
h1=rectangle('Position',[s(1)-W/2,0,W,H],'Curvature',.1,'FaceColor',[1 0.1 0.1],'EdgeColor',[1 1 1]);
    
    % plot the pole
h2=plot([s(1) px(1)],[0 py(1)],'w','LineWidth',2);
h3=rectangle('Position',[px(1)-mr/2,py(1)-mr/2,mr,mr],'Curvature',1,'FaceColor',[.3 0.3 1],'EdgeColor',[1 1 1]);
h4=text(0.95, 1.5, ['phi: ', num2str(1)]);
set(h4,'color','w', 'fontsize', 14);
xlim([-15 5]);
ylim([-5 5]);
set(gca,'Color','k','XColor','w','YColor','w')
set(gcf,'Color','k')

pause(0.1) 
tic

% animation in a for loop 
for k=1:length(t)

  % update pole and cart position
  set(h1, 'position',[s(k)-W/2,0,W,H]);
  set(h2, 'XData',[s(k) px(k)], 'YData', [H/2 py(k)]);
  set(h3, 'position', [px(k)-mr/2,py(k)-mr/2,mr,mr]);
  set(h4, 'string',['time: ', num2str(t(k))]);
 
  drawnow();
  pause(0.1)
  
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
