close all; clear all

imgname = 'trace-without-GP';
outfile = fullfile(pwd,'/simresults/20-01-15-out-GP-without-GP.mat');
kplot = 558;

imgname = 'trace-with-GP';
outfile = fullfile(pwd,'/simresults/20-01-15-out-GP-with-GP-optimized.mat');
kplot = 412

load(outfile)


[trackdata, x0, th0, w] = RaceTrack.loadTrack_02();
track = RaceTrack(trackdata, x0, th0, w);

trackAnim = SingleTrackAnimation(track, out.mu_x_pred_opt, out.var_x_pred_opt, out.u_pred_opt, out.x_ref);
trackAnim.initTrackAnimation();
drawnow;

% k = find(~isnan(out.xhat(1,:)), 1, 'last' ) - 1;

trackAnim.mu_x_pred_opt  = out.mu_x_pred_opt;
trackAnim.var_x_pred_opt = out.var_x_pred_opt;
trackAnim.u_pred_opt     = out.u_pred_opt;
trackAnim.x_ref          = out.x_ref;
trackAnim.updateTrackAnimation(kplot); % 558


delete(trackAnim.h_car)
delete(trackAnim.h_x_ref)
delete(trackAnim.h_mu_x_pred_opt)
cellfun( @(h) delete(h), trackAnim.h_var_x_pred_opt )
legend('off')
title('')
fp.savefig(imgname,'format','epsc')





%% CHECK RELAXED BARRIER FUNCTION PLOT


x = -0.5:0.001:0.5;
% Relaxied barrier function for (x<=lambda)
gamma = 10000;
lambda = -0.2;
y = 0.5*(sqrt((4+gamma*(lambda-x).^2)/gamma) - (lambda-x)); 
figure('Color','w'); hold on; %grid on;
plot(x,y,'LineWidth',2)
ylim([0,0.1])
xlim([-0.4,0])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
%%%%
lambda = -0.2;
x = -0.6:0.001:(lambda-eps);
y = 0.1*-log(lambda-x);
figure('Color','w'); hold on; %grid on;
plot(x,y,'LineWidth',2)
xlim([-0.6,0])

