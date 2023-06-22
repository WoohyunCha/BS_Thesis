clear all; close all ;clc;

%% Simulation Settings
N = 2; % number of bodies
T_end = 1; % simulation end time
T_list = [0.1 0.05 0.02 0.01 0.008 0.004 0.002]; % step sizes to use for TO
loss_list = [];

T = 0.01;

%% Dynamics setting
m = [5;5]; % mass of links
L = [1;1]; % length of links
r = [1;1]; % position of the COM of each link
J = [m(1)*(L(1)^2+0.2^2)/36; m(2)*(L(2)^2+0.2^2)/36]; % Inertia
g = 0;


%% Initialize
time = zeros(1,T_end/T);
M = zeros(2,2*T_end/T); % Mass matrix
C = zeros(2,2*T_end/T); % C&C matrix 
G = zeros(2,T_end/T); % Gravity vector
U = zeros(T_end/T,1); % Torque input vector. [u1]

% initial conditions, q(-1) and q(0)
q_0 = [pi/2;0]; % initial condition q(0)
v_0 = [0;0]; % initial condition q(-1) -> v(0) = 0

%% Simulation setting
simul = false;
test = false; % To check if integrator works properly. Free swing simulation
fig = true;
video = false;
%% TO
if test   
    if ~exist('./test_pmi/', 'dir')
       mkdir('./test_pmi/')
    end
    q_0 = [pi/2;0];
    v_m = [0;0];
    for j = 1:length(T_list) 
        % T_end = 20;
        load(sprintf('./inputs_pmi/U_%f_endtime_%f.mat', T_list(j), T_end), 'U');
        U_simul = U;
        T_simul = T_list(j);
        [qf, M, C, G] = forward(q_0,v_0, r, U_simul, m, L, g, J, T_simul, T_end);
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
    % ylim([1.2, 2.2]);
    title("Trajectory, q1");
    legend('h = 0.1', 'h = 0.05','h = 0.02','h = 0.01','h = 0.008','h = 0.005','h = 0.002')
    saveas(gcf, sprintf('./test_pmi/q1.jpg'));

    figure(2);
    for j = 1:length(T_list)
        t = 0:T_list(j):T_end;
        q = eval(sprintf('q_%d', j));
        plot(t, q(2:2:end));
        hold on;
    end
    % ylim([-1.6, -0.6]);
    title("Trajectory, q2");
    legend('h = 0.1', 'h = 0.05','h = 0.02','h = 0.01','h = 0.008','h = 0.005','h = 0.002')
    saveas(gcf, sprintf('./test_pmi/q2.jpg'));
else
    q_target = [0; 0];
    if simul == true
        if ~exist('./transfer_pmi/', 'dir')
           mkdir('./transfer_pmi/')
        end
        T_simul = 0.002;
        U_simul = zeros(int64(T_end/T_simul),1);
        for j = 1:length(T_list)
            tic
            T = T_list(j);
            load(sprintf('./inputs_pmi/U_%f_endtime_%f.mat', T, T_end), 'U');
            for i = 1:T_end/T
                U_simul(int64(T/T_simul * (i-1)+1) : int64(T/T_simul*i) ) = U(i) ;
            end
            qf = forward(q_0, v_0, r, U_simul, m, L, g, J, T_simul, T_end);
            loss = ((qf(end-1:end,1) - q_target)')*(qf(end-1:end) - q_target) ;
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
            plot(t, qf(2*int64(scale)-1:2*int64(scale):2*int64(T_end/T_simul)-1));
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
            plot(t, qf(2*int64(scale):2*int64(scale):2*int64(T_end/T_simul)));
            hold on;
            plot(t, q_target(2) * ones(length(t)));
            xlabel('time, [s]');
            ylabel('angle, [rad]');
            title('q_2');
            legend('real', 'target');
            saveas(gcf, sprintf('./transfer_pmi/q2_%f_endtime_%f.jpg', T, T_end))
            if (video)
                for k = scale :scale :T_end/T_simul
                    k = int64(k);
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
        save('transfer_pmi/loss_list_pmi.mat', 'loss_list');
        save('transfer_pmi/T_list_pmi.mat', 'T_list');
    else
        if ~exist('./inputs_pmi/', 'dir')
           mkdir('./inputs_pmi/')
        end
        for j = 1: length(T_list)
            tic
            U = zeros(T_end/T_list(j),1); % Torque input vector. [u1]
            q_obj = zeros(2*T_end/T_list(j),1); % Target
            options = optimoptions('fminunc', 'MaxFunctionEvaluations', 200000, 'display', 'off');
            U = fminunc(@(U)ObjFunc(U, q_0, v_0,q_obj,r, m, L, g, J, T_list(j), T_end), U,options);
            U_simul = U;
            T_simul = T_list(j);
            save(sprintf('./inputs_pmi/U_%f_endtime_%f.mat', T_list(j), T_end), 'U');
            [qf, ~, ~, ~] = forward(q_0,v_0, r, U_simul, m, L, g, J, T_simul, T_end);
            toc
            fprintf("TO done for h = %f\n", T_list(j));
                        scale = (T_end/T_simul)/100;
            t = T_list(j):T_list(j):T_end;
            if (fig) 
                figure(3*j-2);
            else
                figure('visible','off');
            end
            plot(t, qf(1:2:end));
            hold on;
            plot(t, q_target(1) * ones(length(t)));
            xlabel('time, [s]');
            ylabel('angle, [rad]');
            title('q_1');
            legend('real', 'target');
            if (fig) 
                figure(3*j-1);
            else
                figure('visible','off');
            end
            plot(t, qf(2:2:end));
            hold on;
            plot(t, q_target(2) * ones(length(t)));
            xlabel('time, [s]');
            ylabel('angle, [rad]');
            title('q_2');
            legend('real', 'target');

        end
    end
end

T = T_simul;
%% Plot


%% Functions
function [q, v, M, C, G] = forward(q_0,v_0, r, U, m, L, g, J, T, T_end)
    % Forward simulation
    q = zeros(2*int64(T_end/T),1);
    v = zeros(2*int64(T_end/T),1);
    for k = 1:int64(T_end/T)
        if (k == 1)
            [q(2*k-1:2*k), v(2*k-1:2*k)] = solve_dynamics(q_0, v_0,r, [U(k);0], m, L, g, J, T);
        else
            [q(2*k-1:2*k), v(2*k-1:2*k)] = solve_dynamics(q(2*k-3:2*k-2), v(2*k-3:2*k-2),r, [U(k);0], m, L, g, J, T);
        end
    % Calculate dynamics matrices
        M(:, 2*k-1:2*k) = Mass(q(2*k-1:2*k),r, m, L, J);
        C(:, 2*k-1:2*k) = Cori(q(2*k-1:2*k), v(2*k-1:2*k) ,r,  m, L);
        G(:, k) = Grav(q(2*k-1:2*k),r, m, L, g);
    end
end

function f = ObjFunc(U, q_0, v_0, q_obj,r, m, L, g, J, T, T_end) % objective function without gradient
    % Calculate objective
    % Forward simulation
    if T <= 0.001
        U_use = zeros(T_end / T, 1);
        for i = 1:T_end/0.002
          U_use((0.002/T * (i-1)+1) : (0.002/T*i) ) =  U(i);
          [q, ~,~,~] = forward(q_0, v_0,r, U_use, m, L, g, J, T, T_end);
        end
    else
        [q, ~,~,~] = forward(q_0, v_0,r, U, m, L, g, J, T, T_end);
    end
    f = norm(q(2*T_end/T-1:2*T_end/T) - q_obj(2*T_end/T-1:2*T_end/T));
end

function [q,v] = solve_dynamics(q_0_, v_0_,r_, u_, m_, L_, g_, J_, T_) % Forward simulation works well
    m = Mass(q_0_,r_, m_, L_, J_);
    sqrt_m = sqrtm(m);
    c = Cori(q_0_, v_0_,r_, m_, L_);
    dm = sylvester(sqrt_m, sqrt_m, c+c');
    Q = sqrt_m\c/sqrt_m - dm/sqrt_m;
    g = Grav(q_0_, r_, m_, L_, g_);
    dg = dg_dq(q_0_, r_,m_, L_, g_);
    z_0_ = sqrt_m * v_0_;
    zeta = (eye(2)/T_ + Q/2 + sqrt_m\dg*T_/4/sqrt_m)\(z_0_/T_ - Q/2*z_0_ + sqrt_m\(u_-g-dg/4*T_/sqrt_m*z_0_));
    q = sqrt_m\(zeta + sqrt_m*v_0_)/2*T_ + q_0_;
    v = sqrtm(Mass(q, r_, m_, L_, J_))\zeta;
end

function v = Grav(q_,r_, m_, L_, g_) % Gravity vector
    v = [m_(1)*g_*r_(1)*L_(1)*sin(q_(1))+m_(2)*g_*L_(1)*sin(q_(1))+m_(2)*g_*r_(2)*L_(2)*sin(q_(1)+q_(2));
        m_(2)*g_*r_(2)*L_(2)*sin(q_(1)+q_(2))];
end

function m = Mass(q_,r_, m_, L_, J_) % mass matrix
    m = [m_(1)*r_(1)^2*L_(1)^2+m_(2)*r_(2)^2*L_(2)^2+m_(2)*L_(1)^2+2*m_(2)*r_(2)*L_(1)*L_(2)*cos(q_(2)), m_(2)*(r_(2)^2*L_(2)^2+r_(2)*L_(1)*L_(2)*cos(q_(2)));
        m_(2)*(r_(2)^2*L_(2)^2+r_(2)*L_(1)*L_(2)*cos(q_(2))), m_(2)*r_(2)^2*L_(2)^2] + [J_(1)+J_(2), J_(2); J_(2), J_(2)];
end

function c = Cori(q_, v_, r_, m_, L_) % C&C matrix
    c= [-m_(2)*r_(2)*L_(1)*L_(2)*sin(q_(2))*v_(2), -m_(2)*r_(2)*L_(1)*L_(2)*sin(q_(2))*v_(1)-m_(2)*r_(2)*L_(1)*L_(2)*sin(q_(2))*v_(2);
        m_(2)*r_(2)*L_(1)*L_(2)*sin(q_(2))*v_(1), 0];
end

function dg = dg_dq(q_,r_, m_, L_, g_) % Jacobian of the gravity vector
    dg = [m_(1)*g_*L_(1)*cos(q_(1))*r_(1)+m_(2)*g_*L_(1)*cos(q_(1))+m_(2)*g_*L_(2)*cos(q_(1)+q_(2))*r_(2), m_(2)*g_*L_(2)*cos(q_(1)+q_(2))*r_(2);
           m_(2)*g_*L_(2)*cos(q_(1)+q_(2))*r_(2), m_(2)*g_*L_(2)*cos(q_(1)+q_(2))*r_(2)];
end
