classdef MPC
    %MPC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        maxiter = 200
        
    end
    
    methods
        function obj = MPC(inputArg1,inputArg2)
            
        end
        
        function outputArg = set(obj,inputArg)
            
        end
    end
end


%% User input

saveimages = true;

if saveimages
    saveimg = @(fig,name,format) fp.savefig(fig,name,format);
else
    saveimg = @(fig,name,format) 0;
end

% cost function parameters
Q = eye(2);
R = 1;

% horizon length
tf = 10;            % horizon length
N = 50;            % number of steps
h = tf/N;           % step size

% continuous dynamic equation parameters
x0 = [-1 1]';
Ac = [1 2; 
      1 1];
Bc = [1;
      2];

n = length(Ac);     % state dimension
m = size(Bc,2);     % input dimension


%% Discretize model

% exact discretization      - c2d(ss(Ac,Bc,eye(n),0), h, 'zoh')
Ad_1 = expm(Ac*h);
syms t
Bd_1 = double( int(expm(Ac*t), t, 0, h) * Bc );


% Euler discretization
Ad_2 = eye(n) + h*Ac;
Bd_2 = Bc*h;


% compare results
Ad_1, Ad_2
Bd_1, Bd_2


%% Question F) Lets verify the difference in the solution when using the two different discretization methods

[H, f, d, Aeq, beq]        = define_opt_control_problem (Ad_1, Bd_1, Q, R, x0, N, n, m);
[x_opt_1, u_opt_1]         = solve_opt_control_quadprog (H, f, d, Aeq, beq, N, n, m);
[x_opt_1_KKT, u_opt_1_KKT] = solve_opt_control_KKT      (H, f, d, Aeq, beq, N, n, m);

norm(x_opt_1-x_opt_1_KKT)

fig = figure('Color','w'); hold on; grid on;
plot_sim ( fig, x0, x_opt_1, u_opt_1, h, tf)
title('Opt. control solution using exact discretization and quadprog solver');
saveimg(fig, 'question-f-exact-quadprog', 'jpg');

fig = figure('Color','w'); hold on; grid on;
plot_sim ( fig, x0, x_opt_1_KKT, u_opt_1_KKT, h, tf)
title('Opt. control solution using exact discretization and KKT-cond');
saveimg(fig, 'question-f-exact-KKT', 'jpg');


[H, f, d, Aeq, beq]        = define_opt_control_problem (Ad_2, Bd_2, Q, R, x0, N, n, m);
[x_opt_2, u_opt_2]         = solve_opt_control_quadprog (H, f, d, Aeq, beq, N, n, m);
[x_opt_2_KKT, u_opt_2_KKT] = solve_opt_control_KKT      (H, f, d, Aeq, beq, N, n, m);

norm(x_opt_2-x_opt_2_KKT)

fig = figure('Color','w'); hold on; grid on;
plot_sim ( fig, x0, x_opt_2, u_opt_2, h, tf)
title('Opt. control solution using Euler discretization and quadprog solver');
saveimg(fig, 'question-f-euler', 'jpg');

fig = figure('Color','w'); hold on; grid on;
plot_sim ( fig, x0, x_opt_2_KKT, u_opt_2_KKT, h, tf)
title('Opt. control solution using Euler discretization and KKT-cond');
saveimg(fig, 'question-f-euler', 'jpg');


fig = figure('Color','w'); hold on; grid on;
plot_sim ( fig, x0*0, x_opt_1-x_opt_2, u_opt_1-u_opt_2, h, tf)
title({'Difference of optimal states and control inputs',
       'for exact and euler discretization methods',
       'using quadprog solver'});
saveimg(fig, 'question-f-error', 'jpg');



%% Question G) What happens if we change Q?

[H, f, d, Aeq, beq] = define_opt_control_problem (Ad_1, Bd_1, 0.2*Q, R, x0, N, n, m);
[x_opt, u_opt]      = solve_opt_control_quadprog (H, f, d, Aeq, beq, N, n, m);

fig = figure('Color','w'); hold on; grid on;
plot_sim ( fig, x0, x_opt, u_opt, h, tf)
title('Opt. control solution using exact discretization - Q = 0.2*I_2');
saveimg(fig, 'question-g-alpha-0p2', 'jpg');


[H, f, d, Aeq, beq] = define_opt_control_problem (Ad_1, Bd_1, 1*Q, R, x0, N, n, m);
[x_opt, u_opt]      = solve_opt_control_quadprog (H, f, d, Aeq, beq, N, n, m);

fig = figure('Color','w'); hold on; grid on;
plot_sim ( fig, x0, x_opt, u_opt, h, tf)
title('Opt. control solution using exact discretization - Q = 1*I_2');
saveimg(fig, 'question-g-alpha-1p0', 'jpg');


[H, f, d, Aeq, beq] = define_opt_control_problem (Ad_1, Bd_1, 40*Q, R, x0, N, n, m);
[x_opt, u_opt]      = solve_opt_control_quadprog (H, f, d, Aeq, beq, N, n, m);

fig = figure('Color','w'); hold on; grid on;
plot_sim ( fig, x0, x_opt, u_opt, h, tf)
title('Opt. control solution using exact discretization - Q = 40*I_2');
saveimg(fig, 'question-g-alpha-40p0', 'jpg');


%% Question H) Compare continuous-time and discrete-time solution using the optimal control input

u     = u_opt_1;
x_opt = x_opt_1;

sys = ss(Ac,Bc,eye(n),0);
[y,t,x] = lsim(sys, u_opt_1, 0:h:tf-h, x0, 'zoh');

figure('Color','w'); hold on; grid on;
plot(1:N-1, x_opt, '-', 'LineWidth',2)
plot(1:N-1, x(2:end,:)', '-.', 'LineWidth',2)

fig = figure('Color','white','Position',[-1734   449   715   282]); 
hold on; grid on;
plot ( 1:N-1, x_opt-x(2:end,:)', 'LineWidth',2)
legend('error x(1)','error x(2)')
title('Error between continuous time and discrete time simulation');
saveimg(fig, 'question-h-discretization-error', 'jpg');


% close all;



%% Helper functions


function [H, f, d, Aeq, beq] = define_opt_control_problem ( Ad, Bd, Q, R, x0, N, n, m)
    %% Define Optimal control matrices

    gamma = zeros(N*n,N*m);
    for i=0:N-2
        gamma = gamma + kron(diag(ones(N-i,1),-i),Ad^i*Bd);     % fill lower-diagonals with Ad^i*Bd
    end
    gamma = gamma(1:(N-1)*n,:);     % remove last constraint that includes x_N
    
    omega = [];
    for i=1:N
        omega = [omega; Ad^i];                              % extend matrix
    end
    omega = omega(1:(N-1)*n,:);     % remove last constraint that includes x_N
    
    Qb = kron(eye(N-1),Q);
    Rb = kron(eye(N),R);                        % Rb = blkdiag(R, R, ..., R)
    
    %% Define quadratic minimization problem

    % J = (xb'*Qb*xb + ub'*Rb*ub + x0'*Q*x0)    =   0.5*[xb,ub]'*H*[xb,ub] + f'*[xb,ub]' + d
    %   where:
    H = 2 * blkdiag(Qb,Rb);
    f = zeros(n*(N-1)+m*N,1);   % f = zeros((n+m)*N,1);
    d = x0' * Q * x0;

    % equality constraints:  
    %  xb = omega*x0 + gamma*ub    =>   Aeq*[xb,ub] = beq
    %   where:
    Aeq = [ eye(n*(N-1)), -gamma ];  % Aeq = [ eye(n*N), -gamma ];
    beq = omega * x0;
end


function [x_opt, u_opt] = solve_opt_control_quadprog (H, f, d, Aeq, beq, N, n, m)
    %% Solve optimal control problem numerically

    % user input
    tol  = 1e-12;
    maxiter = 200;
    
    options = optimoptions('quadprog',...
                           'Algorithm','interior-point-convex',...
                           'ConstraintTolerance',tol,...
                           'MaxIterations', maxiter,...
                           'Display','iter');
    [Y,~,~] = quadprog(H, f, [],[], Aeq,beq, [],[], [],options);
    
    % check constraints error: (Aeq*Y-beq)'

    % reshape solution
    x_opt = reshape( Y(1:(N-1)*n), n, []);
    u_opt = reshape( Y((N-1)*n+1:end), m, []);
end

function [x_opt, u_opt] = solve_opt_control_KKT (H, f, d, Aeq, beq, N, n, m)
    %% Solve optimal control problem using KKT conditions
    
    % 1 step solution
    if rank(H) == size(H,1)
        nu = -(Aeq/H*Aeq')\(Aeq/H*f + beq);
        y  = -H\(Aeq'*nu + f);
        
    % 1 iterative solution - Newton's method
    else    
        % user input
        tol  = 1e-12;
        maxiter = 20;

        y = zeros((n+m)*N,1);
        nu= zeros(n*N,1);

        iter  = 1;
        delta = Inf;
        while iter<=maxiter && norm(delta) >= tol
            A = [ 0.5*(H'+H)  Aeq';
                  Aeq         zeros(n*N,n*N)];
            b = -[0.5*(H'+H)*y+f+Aeq'*nu;
                  Aeq*y-beq];
            delta = A \ b;
            y = y  + delta(1:N*(n+m));
            nu= nu + delta(N*(n+m)+1:end);
            fprintf('iter: %4d, tolerance:%f\n',iter,norm(delta));
            iter = iter+1;
        end
    end
    
    % reshape solution
    x_opt = reshape( y(1:(N-1)*n), n, []);
    u_opt = reshape( y((N-1)*n+1:end), m, []);
end


function plot_sim ( fig, x0, x_opt, u_opt, h, tf )
    %% Show results
    figure(fig);
    clr = lines(20);
    % subplot(3,1,1)
    stairs(0:h:tf-h, [x0(1), x_opt(1,:)]', '-','Color',clr(1,:),'LineWidth',2,'DisplayName','x(1)');
    stairs(0:h:tf-h, [x0(2), x_opt(2,:)]', '-','Color',clr(2,:),'LineWidth',2,'DisplayName','x(2)');
    stairs(0:h:tf-h, u_opt(1,:)',          '-','Color',clr(3,:),'LineWidth',2,'DisplayName','u'); 

    legend
    xlabel('time [s]')
end

