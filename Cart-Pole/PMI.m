clear all; close all ;clc;

%% Simulation Settings
N = 2; % number of bodies
T_end = 1; % simulation end time
T_list = [0.05 0.02 0.01 0.008 0.004 0.002];
loss_list = [];

%% Dynamics setting
T_test = 0.01;
m = [5;1]; % mass of links
L = 0.5; % length of links
J = 1/2*m(2)*0.1^2; % Inertia
g = 9.8;


%% Initialize

% initial conditions, q(-1) and q(0)
q_0 = [0;0]; % initial condition q(0)
v_0 = [0;0]; % initial condition q(-1) -> v(0) = 0

%% Simulation setting
simul = false;
test = false; % To check if integrator works properly. Free swing simulation
video = false;
fig = false;
%% TO
if test
    if ~exist('./test_pmi/', 'dir')
       mkdir('./test_pmi/')
    end
    q_0 = [0;0];
    v_m = [0;0];
    for j = 1:length(T_list) 
        % T_end = 20;
        load(sprintf('./inputs_pmi/U_%f_endtime_%f.mat', T_list(j), T_end), 'U');
        U_simul = U;
        T_simul = T_list(j);
        [qf, M, C, G] = forward(q_0,v_0, U_simul, m, L, g, J, T_simul, T_end);
        fieldValue = [q_0;qf];
        assignin('base', sprintf('q_%d', j), fieldValue);
        save(sprintf('./test_pmi/q_%f_endtime_%f.mat', T_list(j), T_end), sprintf('q_%d', j));
    end
    figure(1);
    for j = 1:length(T_list)
        t = 0:T_list(j):T_end;
        q = eval(sprintf('q_%d', j));
        plot(t, q(1:2:end));
        hold on;
    end
    title("Trajectory, q1");
    legend('h = 0.05', 'h = 0.02', 'h = 0.01', 'h = 0.008', 'h = 0.004','h = 0.002')
    saveas(gcf, sprintf('./test_pmi/q1.jpg'));

    figure(2);
    for j = 1:length(T_list)
        t = 0:T_list(j):T_end;
        q = eval(sprintf('q_%d', j));
        plot(t, q(2:2:end));
        hold on;
    end
    ylim([-1.5,0]);
    title("Trajectory, q2");
    legend('h = 0.05', 'h = 0.02', 'h = 0.01', 'h = 0.008', 'h = 0.004','h = 0.002')
    saveas(gcf, sprintf('./test_pmi/q2.jpg'));
else
    q_target = [0; pi];
    if simul == true
        if ~exist('./transfer_pmi/', 'dir')
           mkdir('./transfer_pmi/')
        end
        T_simul = 0.0001;
        U_simul = zeros(T_end/T_simul,1);
        for j = 1:length(T_list)
            tic
            T = T_list(j);
            load(sprintf('./inputs_pmi/U_%f_endtime_%f.mat', T, T_end), 'U');
            for i = 1:T_end/T
                U_simul((T/T_simul * (i-1)+1) : (T/T_simul*i) ) = U(i) ;%* T_simul/T;
            end
            qf = forward(q_0, v_0, U_simul, m, L, g, J, T_simul, T_end);
            loss = ((qf(2*T_end/T_simul-1:2*T_end/T_simul,1) - q_target)')*(qf(2*T_end/T_simul-1:2*T_end/T_simul) - q_target) ;
            save(sprintf('./transfer_pmi/q_%f_siumultime_%f_endtime_%f.mat', T, T_simul, T_end), 'qf');
            save(sprintf('./transfer_pmi/loss_%f_simultime_%f_endtime_%f.mat', T, T_simul, T_end), 'loss');
            loss_list = [loss_list, loss];
            toc
            fprintf("Simulation done for h = %f\n", T_list(j));
            
            scale = (T_end/T_simul)/100;
            t = T_simul*scale:T_simul*scale:T_end;
            if (fig) 
                figure(3*j-2);
            else
                figure('visible','off');
            end
            plot(t, qf(2*scale-1:2*scale:2*T_end/T_simul-1));
            hold on;
            plot(t, q_target(1) * ones(length(t)));
            xlabel('time, [s]');
            ylabel('angle, [rad]');
            title('q_1');
            legend('real', 'target');
            saveas(gcf, sprintf('./transfer_pmi/q1_%f_endtime_%f.jpg', T, T_end));
            if (fig) 
                figure(3*j-1);
            else
                figure('visible','off');
            end
            plot(t, qf(2*scale:2*scale:2*T_end/T_simul));
            hold on;
            plot(t, q_target(2) * ones(length(t)));
            xlabel('time, [s]');
            ylabel('angle, [rad]');
            title('q_2');
            legend('real', 'target');
            saveas(gcf, sprintf('./transfer_pmi/q2_%f_endtime_%f.jpg', T, T_end))
            if (video)
                for k = scale :scale :T_end/T_simul
                    figure(3*j)
                    q1 = qf(2*k-1);
                    q2 = qf(2*k);
                    x = zeros(2, 2);
                    x(:, 1) = [q1 ; 0];
                    x(:, 2) = x(:, 1) + L*[sin(q2); -cos(q2)];
                
                    clf;
                    plot(x(1, :), x(2, :), 'o-');
                    hold on;
                    plot(x(1, 1), x(2, 1), 'square');
                    hold on;
                    plot([-10,10], [0,0], '-');
                    axis equal;
                    xlim([-2, 2]);
                    ylim([-1,1]);
                    drawnow;
                    frame = getframe(gcf);
                end
            end
        end
    else
        if ~exist('./inputs_pmi/', 'dir')
           mkdir('./inputs_pmi/')
        end
        for j = 1: length(T_list)
            tic
            U = zeros(T_end/T_list(j),1); % Torque input vector. [u1]
            options = optimoptions('fminunc', 'MaxFunctionEvaluations', 200000, 'display', 'off');
            U = fminunc(@(U)ObjFunc(U, q_0, v_0,q_target, m, L, g, J, T_list(j), T_end), U,options);
            U_simul = U;
            T_simul = T_list(j);
            save(sprintf('./inputs_pmi/U_%f_endtime_%f.mat', T_list(j), T_end), 'U');
            [qf, ~, ~, ~] = forward(q_0,v_0, U_simul, m, L, g, J, T_simul, T_end);
            save(sprintf('./inputs_pmi/q_%f_endtime_%f.mat', T_list(j), T_end), 'qf');
            toc
            fprintf("TO done for h = %f\n", T_list(j));
            
            scale = (T_end/T_simul)/10;
            t = T_simul*scale:T_simul*scale:T_end;
            if (fig) 
                figure(3*j-2);
            else
                figure('visible','off');
            end
            plot(t, qf(2*scale-1:2*scale:2*T_end/T_simul-1));
            hold on;
            plot(t, q_target(1) * ones(length(t)));
            xlabel('time, [s]');
            ylabel('angle, [rad]');
            title('q_1');
            legend('real', 'target');
            saveas(gcf, sprintf('./inputs_pmi/q1_%f_endtime_%f.jpg', T_simul, T_end));
            if (fig) 
                figure(3*j-1);
            else
                figure('visible','off');
            end
            plot(t, qf(2*scale:2*scale:2*T_end/T_simul));
            hold on;
            plot(t, q_target(2) * ones(length(t)));
            xlabel('time, [s]');
            ylabel('angle, [rad]');
            title('q_2');
            legend('real', 'target');
            saveas(gcf, sprintf('./inputs_pmi/q2_%f_endtime_%f.jpg', T_simul, T_end))
            if (video)
                for k = scale : scale :T_end/T_simul
                    figure(3*j)
                    q1 = qf(2*k-1);
                    q2 = qf(2*k);
                    x = zeros(2, 2);
                    x(:, 1) = [q1 ; 0];
                    x(:, 2) = x(:, 1) + L*[sin(q2); -cos(q2)];
                
                    clf;
                    plot(x(1, :), x(2, :), 'o-');
                    hold on;
                    plot(x(1, 1), x(2, 1), 'square');
                    hold on;
                    plot([-10,10], [0,0], '-');
                    axis equal;
                    xlim([-2, 2]);
                    ylim([-1,1]);
                    drawnow;
                    frame = getframe(gcf);
                end
            end
        end
    end
end

T = T_simul;
%% Plot




%% Functions
function [q, M, C, G] = forward(q_0,v_0, U, m, L, g, J, T, T_end)
    % Forward simulation
    q = zeros(2*T_end/T,1);
    v = zeros(2*T_end/T,1);
    for k = 1:T_end/T
        if (k == 1)
            [q(2*k-1:2*k), v(2*k-1:2*k)] = solve_dynamics(q_0, v_0, [U(k);0], m, L, g, J, T);
        else
            [q(2*k-1:2*k), v(2*k-1:2*k)] = solve_dynamics(q(2*k-3:2*k-2), v(2*k-3:2*k-2), [U(k);0], m, L, g, J, T);
        end
    % Calculate dynamics matrices
        M(:, 2*k-1:2*k) = Mass(q(2*k-1:2*k),m, L, J);
        C(:, 2*k-1:2*k) = Cori(q(2*k-1:2*k), v(2*k-1:2*k) , m, L);
        G(:, k) = Grav(q(2*k-1:2*k), m, L, g);
    end
end

function f = ObjFunc(U, q_0, q_m, q_target, m, L, g, J, T, T_end) % objective function without gradient
    % Calculate objective
    % Forward simulation
    if T <= 0.001
        U_use = zeros(T_end / T, 1);
        for i = 1:T_end/0.002
          U_use((0.002/T * (i-1)+1) : (0.002/T*i) ) =  U(i);
          [q, ~,~,~] = forward(q_0, q_m, U_use, m, L, g, J, T, T_end);
        end
    else
        [q,~,~,~] = forward(q_0, q_m, U, m, L, g, J, T, T_end);
    end
    f = ((q(2*T_end/T-1:2*T_end/T,1) - q_target)')*(q(2*T_end/T-1:2*T_end/T) - q_target) + (1e-5)*norm(U)^2;
end

function [q,v] = solve_dynamics(q_0_, v_0_, u_, m_, L_, g_, J_, T_) % Forward simulation works well
    m = Mass(q_0_, m_, L_, J_);
    sqrt_m = sqrtm(m);
    c = Cori(q_0_, v_0_, m_, L_);
    dm = sylvester(sqrt_m, sqrt_m, c+c');
    Q = sqrt_m\c/sqrt_m - dm/sqrt_m;
    g = Grav(q_0_, m_, L_, g_);
    dg = dg_dq(q_0_,m_, L_, g_);
    z_0_ = sqrt_m * v_0_;
    zeta = (eye(2)/T_ + Q/2 + sqrt_m\dg*T_/4/sqrt_m)\(z_0_/T_ - Q/2*z_0_ + sqrt_m\(u_-g-dg/4*T_/sqrt_m*z_0_));
    q = sqrt_m\(zeta + sqrt_m*v_0_)/2*T_ + q_0_;
    v = sqrtm(Mass(q, m_, L_, J_))\zeta;
end

function v = Grav(q_, m_, L_, g_) % Gravity vector
    v = [0; 
        m_(2)*g_*L_*sin(q_(2))];
end

function m = Mass(q_, m_, L_, J_) % mass matrix
    m = [m_(1)+m_(2), m_(2) * L_*cos(q_(2));
        m_(2)*L_*cos(q_(2)), m_(2)*L_^2+J_];
end

function c = Cori(q_, v_, m_, L_)
    c = [0 -m_(2)*L_*sin(q_(2))*v_(2);
        0 0];
end

function dg = dg_dq(q_, m_, L_, g_) % Jacobian of the gravity vector
    dg = [0, 0;
        0, m_(2)*g_*L_*cos(q_(2))];
end

