
clear all; close all; clc;

% Linear wheel dynamics from nominal model
c_f = MotionModelGP_SingleTrack_nominal.c_f;
c_r = MotionModelGP_SingleTrack_nominal.c_r;


% Pacejka lateral dynamics parameters
B_f = MotionModelGP_SingleTrack_true.B_f;             % stiffnes factor (Pacejka) (front wheel)
C_f = MotionModelGP_SingleTrack_true.C_f;             % shape factor (Pacejka) (front wheel)
D_f = MotionModelGP_SingleTrack_true.D_f;             % peak value (Pacejka) (front wheel)
E_f = MotionModelGP_SingleTrack_true.E_f;             % curvature factor (Pacejka) (front wheel)

B_r = MotionModelGP_SingleTrack_true.B_r;             % stiffnes factor (Pacejka) (rear wheel)
C_r = MotionModelGP_SingleTrack_true.C_r;             % shape factor (Pacejka) (rear wheel)
D_r = MotionModelGP_SingleTrack_true.D_r;             % peak value (Pacejka) (rear wheel)
E_r = MotionModelGP_SingleTrack_true.E_r;             % curvature factor (Pacejka) (rear wheel)

a_r = deg2rad(-30:0.1:30);
a_f = deg2rad(-30:0.1:30);
W_Fy_r = D_r*sin(C_r*atan(B_r*a_r-E_r*(B_r*a_r -atan(B_r*a_r)))); % rear lateral force
W_Fy_f = D_f*sin(C_f*atan(B_f*a_f-E_f*(B_f*a_f -atan(B_f*a_f)))); % front lateral force

figure('Color','w'); hold on; grid on;
plot(rad2deg(a_r),W_Fy_r/1000,'DisplayName','Pacejka tyre model')
plot(rad2deg(a_r),a_r*c_r/1000,'DisplayName','Constant cornering stiffness model')
% title('Rear tyre')
xlabel('Slip angle [deg]');
ylabel('Tyre lateral force [kN]')
legend('Location','northwest')
xlim([-30,30])
ylim([-6,6])
fp.savefig('rear_tyre')

figure('Color','w'); hold on; grid on;
plot(rad2deg(a_f),W_Fy_f/1000,'DisplayName','Pacejka tyre model')
plot(rad2deg(a_f),a_r*c_f/1000,'DisplayName','Constant cornering stiffness model')
% title('Front tyre')
xlabel('Slip angle [deg]');
ylabel('Tyre lateral force [kN]')
legend('Location','northwest')
xlim([-30,30])
ylim([-6,6])
fp.savefig('front_tyre')