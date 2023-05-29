clear all; close all ;clc;

%% Simulation Settings
N = 2; % number of bodies
T_end = 1; % simulation end time
T_list = [0.2 0.1 0.05 0.02 0.01 0.008 0.005 0.002]; % time intervals
loss_list = [];

T = 0.01;

%% Dynamics setting
m = [5;1]; % mass of links
L = 0.5; % length of links
J = 1/2*m(2)*0.1^2; % Inertia
g = 9.8;


%% Initialize
time = zeros(1,T_end/T);
M = zeros(2,2*T_end/T); % Mass matrix
C = zeros(2,2*T_end/T); % C&C matrix 
G = zeros(2,T_end/T); % Gravity vector
U = zeros(T_end/T,1); % Torque input vector. [u1]

% initial conditions, q(-1) and q(0)
q_0 = [0;0]; % initial condition q(0)
v_0 = [0;0]; % initial condition q(-1) -> v(0) = 0

%% Simulation setting
simul = true;
test = false; % To check if integrator works properly. Free swing simulation
%% TO
if test
    U_simul = zeros(2*T_end/T, 1);
    T_simul = T;
    q_obj = zeros(2*T_end/T_simul,1); % Targe
    [qf, vf, ~, ~, ~] = forward(q_0,v_0, r, U_simul, m, L, g, J, T_simul, T_end);
else
    q_target = [0; pi/2];
    if simul == true
        if ~exist('./transfer_pmi/', 'dir')
           mkdir('./transfer_pmi/')
        end
        T_simul = 0.0001;
        U_simul = zeros(T_end/T_simul,1);
        for j = 1:length(T_list)
            T = T_list(j);
            load(sprintf('./inputs_pmi/U_%f_endtime_%f.mat', T, T_end), 'U');
            for i = 1:T_end/T
                U_simul((T/T_simul * (i-1)+1) : (T/T_simul*i) ) = U(i) ;%* T_simul/T;
            end
            loss = ObjFunc(U_simul, q_0, q_m, q_target, m, L, g, J, T_simul, T_end);
            save(sprintf('./transfer_pmi/q_%f_siumultime_%f_endtime_%f.mat', T, T_simul, T_end), 'qf');
            save(sprintf('./transfer_pmi/loss_%f_simultime_%f_endtime_%f.mat', T, T_simul, T_end), 'loss');
            loss_list = [loss_list, loss];
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
            toc
            fprintf("TO done for h = %f\n", T_list(j));
        end
    end
end

T = T_simul;
%% Plot

if (test)
    t = T:T:T_end;
    plot(t, qf(1:2:2*T_end/T-1));
    hold on;
    plot(t, q_obj(1:2:2*T_end/T-1));
    xlabel('time, [s]');
    ylabel('angle, [rad]');
    title('q_1');
    legend('real', 'target');
    figure(2)
    plot(t, qf(2:2:2*T_end/T));
    hold on;
    plot(t, q_obj(2:2:2*T_end/T));
    xlabel('time, [s]');
    ylabel('angle, [rad]');
    title('q_2');
    legend('real', 'target');
    for k = T_end/T /100:T_end/T /100 :T_end/T
        figure(3);
        q1 = qf(2*k-1);
        q2 = qf(2*k);
        x = zeros(2, 2);
        x(:, 1) = [q1 ; 0];
        x(:, 2) = x(:, 1) + L*[sin(q2); -cos(q2)];
    
        clf;
        plot(x(1, :), x(2, :), 'o-');
        axis equal;
        xlim([-1, 1]);
        ylim([-2.5,2.5]);
        drawnow;
        frame = getframe(gcf);
    end
elseif (simul)
        figure(4)
        semilogx(T_list, loss_list);
        xlabel("time step size, [s]");
        ylabel("loss");
        title('PMI integration');
        save('./transfer_PMI/loss_list_PMI', 'loss_list');
        save('./transfer_PMI/T_list_PMI', 'T_list');
end
 


%% Functions
function [q, M, C, G] = forward(q_0,q_m, U, m, L, g, J, T, T_end)
            % Forward simulation
    q = zeros(2*T_end/T,1);
    for k = 1:T_end/T
        if k == 1
            q(2*k-1:2*k) = solve_dynamics(q_0, q_m, [U(k);0], m, L, g, J, T);
        elseif k == 2
            q(2*k-1:2*k) = solve_dynamics(q(2*k-3:2*k-2), q_0, [U(k);0], m, L, g, J, T);
        else
            q(2*k-1:2*k) = solve_dynamics(q(2*k-3:2*k-2), q(2*k-5:2*k-4), [U(k);0], m, L, g,J, T);
        end
    % Calculate dynamics matrices
        M(:, 2*k-1:2*k) = Mass(q(2*k-1:2*k), m, L, J);
        if k == 1
            C(:, 2*k-1:2*k) = Cori(q(2*k-1:2*k), (q(2*k-1:2*k)-q_0)/T,  m, L);
        else
            C(:, 2*k-1:2*k) = Cori(q(2*k-1:2*k), (q(2*k-1:2*k) - q(2*k-3:2*k-2))/T,  m, L);
        end 
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
    f = ((q(2*T_end/T-1:2*T_end/T,1) - q_target)')*(q(2*T_end/T-1:2*T_end/T) - q_target) ;
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

