%% plot tire dynamics

%Pacejka
            B_f = 0.4;              % stiffnes factor (Pacejka) (front wheel)
            C_f = 8;                % shape factor (Pacejka) (front wheel)
            D_f = 4560.4;           % peak value (Pacejka) (front wheel)
            E_f = -0.5;             % curvature factor (Pacejka) (front wheel)
            
            B_r = 0.45;             % stiffnes factor (Pacejka) (rear wheel)
            C_r = 8;                % shape factor (Pacejka) (rear wheel)
            D_r = 4000;             % peak value (Pacejka) (rear wheel)
            E_r = -0.5;             % curvature factor (Pacejka) (rear wheel)
            
% linear
        c_f = 14000 * 2.5 % = 1*g*M/deltamax  % front coornering stiffness (C*delta=Fy~M*a)
        c_r = 14000 * 2.5 % = 2*g*M/deltamax  % rear coornering stiffness
            
            a_r = deg2rad(-25:0.1:25);
            a_f = deg2rad(-25:0.1:25);
            W_Fy_r = D_r*sin(C_r*atan(B_r*a_r-E_r*(B_r*a_r -atan(B_r*a_r)))); % rear lateral force
            W_Fy_f = D_f*sin(C_f*atan(B_f*a_f-E_f*(B_f*a_f -atan(B_f*a_f)))); % front lateral force
            
            figure('Color','w'); hold on; grid on;
            plot(rad2deg(a_r),W_Fy_r/1000,'DisplayName','Pacejka tyre model')
            plot(rad2deg(a_r),a_r*c_r/1000,'DisplayName','Constant coornering stiffness model')
            title('Rear tyre')
            xlabel('Slip angle [deg]');
            ylabel('Tyre lateral force [kN]')
            legend

            figure; hold on; grid on;
            plot(rad2deg(a_f),W_Fy_f/1000,'DisplayName','Pacejka tyre model')
            plot(rad2deg(a_f),a_r*c_f/1000,'DisplayName','Constant coornering stiffness model')
            title('Front tyre')
            xlabel('Slip angle [deg]');
            ylabel('Tyre lateral force [kN]')
            legend