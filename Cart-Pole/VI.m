clear all; close all; clc;

%%% WHEN USING SYMS, MOVE THESE LINES UNDER DYNAMICS SETTING
% syms L(q1x, q1y, q0x, q0y, T_);
% L(q1x, q1y, q0x, q0y, T_) = 1/2*[(q1x-q0x)/T_; (q1y-q0y)/T_]'*Mass([(q1x+q0x)/2; (q1y+q0y)/2], m_, L_, J_)*[(q1x-q0x)/T_; (q1y-q0y)/T_]...
%                               -m_(2)*g_*L_*(1-cos((q1y+q0y)/2));
% 
% dminus = [diff(L, q0x), diff(L, q0y)];
% dplus = [diff(L, q1x), diff(L, q1y)];


%% Simulation Settings
N = 2; % number of bodies
T_end = 1; % simulation end time

T_list = [0.1 0.05 0.02 0.01 0.008 0.005 0.002]; % step sizes to use for TO
loss_list = [];

%% Dynamics setting
m = [5;1]; % mass of links
L = 0.5; % length of links
J = 1/2*m(2)*0.1^2; % Inertia
g = 9.8;

% Variables for symbolic calculations
m_ = m; 
L_ = L;
J_ = J; 
g_ = g;

T_test = 0.01;

% Set modes
test = false;
simul = false; % false means solving TO to get optimal inputs under step sizes of T_list
video = false;
fig = false;
%% Initialize
% initial conditions, q(-1) and q(0)
q_0 = [0;0]; % initial condition q(0)
q_m = [0;0]; % initial condition q(-1) -> v(0) = 0

%% TO
if test
    U_simul = zeros(2*T_end/T_test, 1);
    T_simul = T_test;
    q_obj = zeros(2*T_end/T_simul,1); % Targe
    qf = forward2(q_0, q_m, U_simul, T_simul, T_end);
else
    q_target = [0; pi];
    if simul == true
        if ~exist('./transfer_vi/', 'dir')
           mkdir('./transfer_vi/')
        end
        T_simul = 0.0001;
        U_simul = zeros(T_end/T_simul,1);
        for j = 1:length(T_list)
            tic
            T = T_list(j);
            load(sprintf('./inputs_vi/U_%f_endtime_%f.mat', T, T_end), 'U');
            for i = 1:T_end/T
                U_simul((T/T_simul * (i-1)+1) : (T/T_simul*i) ) = U(i) ;%* T_simul/T;
            end
            qf = forward2(q_0, q_m, U_simul, T_simul, T_end);
            loss = ((qf(2*T_end/T_simul-1:2*T_end/T_simul,1) - q_target)')*(qf(2*T_end/T_simul-1:2*T_end/T_simul) - q_target) ;
            save(sprintf('./transfer_vi/q_%f_siumultime_%f_endtime_%f.mat', T, T_simul, T_end), 'qf');
            save(sprintf('./transfer_vi/loss_%f_simultime_%f_endtime_%f.mat', T, T_simul, T_end), 'loss');
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
            saveas(gcf, sprintf('./transfer_vi/q1_%f_endtime_%f.jpg', T, T_end));
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
            saveas(gcf, sprintf('./transfer_vi/q2_%f_endtime_%f.jpg', T, T_end))
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
        if ~exist('./inputs_vi/', 'dir')
           mkdir('./inputs_vi/')
        end
        for j = 1: length(T_list)
            tic
            U = zeros(T_end/T_list(j),1); % Torque input vector. [u1]
            options = optimoptions('fminunc', 'MaxFunctionEvaluations', 200000, 'display', 'off');
            U = fminunc(@(U)ObjFunc2(U, q_0, q_m,q_target, T_list(j), T_end), U,options);
            U_simul = U;
            T_simul = T_list(j);
            save(sprintf('./inputs_vi/U_%f_endtime_%f.mat', T_list(j), T_end), 'U');
            qf = forward2(q_0,q_m, U_simul, T_simul, T_end);
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
            saveas(gcf, sprintf('./inputs_vi/q1_%f_endtime_%f.jpg', T_simul, T_end));
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
            saveas(gcf, sprintf('./inputs_vi/q2_%f_endtime_%f.jpg', T_simul, T_end))
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

if test
    t = T_test:T_test:T_end;
    plot(t, qf(1:2:2*T_end/T_test-1));
    hold on;
    plot(t, q_obj(1:2:2*T_end/T_test-1));
    xlabel('time, [s]');
    ylabel('angle, [rad]');
    title('q_1');
    legend('real', 'target');
    figure(2)
    plot(t, qf(2:2:2*T_end/T_test));
    hold on;
    plot(t, q_obj(2:2:2*T_end/T_test));
    xlabel('time, [s]');
    ylabel('angle, [rad]');
    title('q_2');
    legend('real', 'target');
    
    v = VideoWriter(sprintf('VI_%f_simultime_%f.avi', T_test, T_simul));
    open(v);
    for k = T_end/T_test /100:T_end/T_test /100 :T_end/T_test
        figure(3)
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
       writeVideo(v,frame);
    end
    close(v);
elseif (simul)
        figure(3*length(T_list)+1)
        semilogx(T_list, loss_list);
        xlabel("time step size, [s]");
        ylabel("loss");
        title('vi integration');
        save('./transfer_vi/loss_list_vi', 'loss_list');
        save('./transfer_vi/T_list_vi', 'T_list');
end

%% Functions

function q = forward2(q_0,q_m, U, T, T_end)
% Forward simulation
    for k = 1:T_end/T
        if k == 1
            q(2*k-1:2*k, 1) = solve_dynamics2(q_0, q_m, [0;0], [U(k);0], T);
        elseif k == 2
            q(2*k-1:2*k, 1) = solve_dynamics2(q(2*k-3:2*k-2), q_0, [U(k-1);0], [U(k); 0], T);
        else
            q(2*k-1:2*k, 1) = solve_dynamics2(q(2*k-3:2*k-2), q(2*k-5:2*k-4), [U(k-1);0], [U(k); 0], T);
        end
    end
end

function f = ObjFunc2(U, q_0, q_m, q_target, T, T_end) % objective function without gradient
    % Calculate objective
    % Forward simulation
    if T <= 0.001
        U_use = zeros(T_end / T, 1);
        for i = 1:T_end/0.002
          U_use((0.002/T * (i-1)+1) : (0.002/T*i) ) =  U(i);
          q = forward2(q_0, q_m, U_use,T, T_end);
        end
    else
        q = forward2(q_0, q_m, U, T, T_end);
    end
    f = ((q(2*T_end/T-1:2*T_end/T,1) - q_target)')*(q(2*T_end/T-1:2*T_end/T) - q_target) ;
end

function q = solve_dynamics2(q_, q_0, u_plus, u_minus, T_)
    q = fsolve(@(x)varint2(x, q_, q_0, u_plus, u_minus, T_), q_, optimoptions('fsolve', 'display', 'off'));
  
end
function f = varint2(x, q_, q_0, u_plus, u_minus, T_) 
q1x = x(1);
q1y = x(2);
q0x = q_(1);
q0y = q_(2);
dm = [((3*(conj(q0x) - conj(q1x)))/conj(T_) + (cos(q0y/2 + q1y/2)*(conj(q0y) - conj(q1y)))/(4*conj(T_)))/T_ + (3*(q0x - q1x))/(T_*conj(T_)) + (cos(q0y/2 + q1y/2)*(q0y - q1y))/(4*T_*conj(T_)), ((51*(conj(q0y) - conj(q1y)))/(400*conj(T_)) + (cos(q0y/2 + q1y/2)*(conj(q0x) - conj(q1x)))/(4*conj(T_)))/T_ - (49*sin(q0y/2 + q1y/2))/20 + ((cos(q0y/2 + q1y/2)/(4*conj(T_)) - (sin(q0y/2 + q1y/2)*(conj(q0y) - conj(q1y)))/(8*conj(T_)))*(q0x - q1x))/T_ + ((51/(400*conj(T_)) - (sin(q0y/2 + q1y/2)*(conj(q0x) - conj(q1x)))/(8*conj(T_)))*(q0y - q1y))/T_];
q1x = q_(1);
q1y = q_(2);
q0x = q_0(1);
q0y = q_0(2);
dp = [- ((3*(conj(q0x) - conj(q1x)))/conj(T_) + (cos(q0y/2 + q1y/2)*(conj(q0y) - conj(q1y)))/(4*conj(T_)))/T_ - (3*(q0x - q1x))/(T_*conj(T_)) - (cos(q0y/2 + q1y/2)*(q0y - q1y))/(4*T_*conj(T_)), - (49*sin(q0y/2 + q1y/2))/20 - ((51*(conj(q0y) - conj(q1y)))/(400*conj(T_)) + (cos(q0y/2 + q1y/2)*(conj(q0x) - conj(q1x)))/(4*conj(T_)))/T_ - ((cos(q0y/2 + q1y/2)/(4*conj(T_)) + (sin(q0y/2 + q1y/2)*(conj(q0y) - conj(q1y)))/(8*conj(T_)))*(q0x - q1x))/T_ - ((51/(400*conj(T_)) + (sin(q0y/2 + q1y/2)*(conj(q0x) - conj(q1x)))/(8*conj(T_)))*(q0y - q1y))/T_];
f = double(T_*(dp + dm) + (u_plus+u_minus)');
end

function m = Mass(q_, m_, L_, J_) % mass matrix
    m = [m_(1)+m_(2), m_(2) * L_*cos(q_(2));
        m_(2)*L_*cos(q_(2)), m_(2)*L_^2+J_];
end
