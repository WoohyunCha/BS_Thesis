clear all; close all; clc;

%%% WHEN USING SYMS, MOVE THESE LINES UNDER DYNAMICS SETTING
% syms L(q1x, q1y, q0x, q0y, T_);
% L(q1x, q1y, q0x, q0y, T_) = 1/2*[(q1x-q0x)/T_; (q1y-q0y)/T_]'*Mass([(q1x+q0x)/2; (q1y+q0y)/2], r_, m_, L_, J_)*[(q1x-q0x)/T_; (q1y-q0y)/T_]...
%                         +(m_(1)*g_*r(1)*L_(1)*cos((q1x+q0x)/2) + m_(2)*g_*(L_(1)*cos((q1x+q0x)/2)+r_(2)*L_(2)*cos((q1x+q0x+q1y+q0y)/2)));
% 
% dminus = [diff(L, q0x), diff(L, q0y)];
% dplus = [diff(L, q1x), diff(L, q1y)];


%% Simulation Settings
N = 2; % number of bodies
T_end = 1; % simulation end time

T_list = [0.1 0.05 0.02 0.01 0.008 0.004 0.002]; % step sizes to use for TO
loss_list = [];

%% Dynamics setting
m = [5;5]; % mass of links
L = [1;1]; % length of links
r = [1;1]; % position of the COM of each link
J = [m(1)*(L(1)^2+0.2^2)/36; m(2)*(L(2)^2+0.2^2)/36]; % Inertia
g = 9.8;

% Variables for symbolic calculations
m_ = m; 
L_ = L; 
r_ = r; 
J_ = J; 
g_ = g;


simul = true; % false means solving variational TO to get optimal inputs under step sizes of T_list
                % True means input transfer
test = false; % To test if the integrator works
fig = false; % To visualize plot
video = false; % To visulaize video
%% Initialize
% initial conditions, q(-1) and q(0)
q_0 = [pi/2;0]; % initial condition q(0)
q_m = [pi/2;0]; % initial condition q(-1) -> v(0) = 0

%% TO
if test
    if ~exist('./test_variational/', 'dir')
       mkdir('./test_variational/')
    end
    q_0 = [pi/2;0];
    q_m = [pi/2;0];
    for j = 1:length(T_list) 
        load(sprintf('./inputs_variational/U_%f.mat', T_list(j)), 'U');
        U_simul = U;
        T_simul = T_list(j);
        qf = forward2(q_0,q_m, U_simul, T_simul, T_end);
        fieldValue = [q_0;qf];
        assignin('base', sprintf('q_%d', j), fieldValue);
        save(sprintf('./test_variational/q_%f_endtime_%f.mat', T_list(j), T_end), sprintf('q_%d', j));
    end
    figure(1);
    for j = 1:length(T_list)
        t = 0:T_list(j):T_end;
        q = eval(sprintf('q_%d', j));
        plot(t, q(1:2:end));
        hold on;
    end
    title("Trajectory, q1");
    legend('h = 0.1', 'h = 0.05','h = 0.02','h = 0.01','h = 0.008','h = 0.005','h = 0.002')
    saveas(gcf, sprintf('./test_variational/q1.jpg'));

    figure(2);
    for j = 1:length(T_list)
        t = 0:T_list(j):T_end;
        q = eval(sprintf('q_%d', j));
        plot(t, q(2:2:end));
        hold on;
    end
    title("Trajectory, q2");
    legend('h = 0.1', 'h = 0.05','h = 0.02','h = 0.01','h = 0.008','h = 0.005','h = 0.002')
    saveas(gcf, sprintf('./test_variational/q2.jpg'));
else
    q_target = [0;0];
    if simul == true % Input transfer
        if ~exist('./transfer_variational/', 'dir')
           mkdir('./transfer_variational/')
        end
        T_simul = 0.00001; % Time step size of 'real world' simulation model
        U_simul = zeros(int64(T_end/T_simul),1);
        for j = 1:length(T_list)
            tic
            T = T_list(j);
            load(sprintf('./inputs_variational/U_%f.mat', T), 'U');
            for i = 1:T_end/T
                U_simul(int64(T/T_simul * (i-1)+1) : int64(T/T_simul*i) ) = U(i) * T_simul/T;
            end
            qf = forward2(q_0, q_m, U_simul, T_simul, T_end);
            loss = ((qf(end-1:end) - q_target)')*(qf(end-1:end) - q_target) ;
            save(sprintf('./transfer_variational/q_%f_siumultime_%f_endtime_%f.mat', T, T_simul, T_end), 'qf');
            save(sprintf('./transfer_variational/loss_%f_simultime_%f_endtime_%f.mat', T, T_simul, T_end), 'loss');
            loss_list = [loss_list, loss];
            toc
            fprintf("Simulation done for h = %f\n", T_list(j));
            %%%PLOT%%%
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
            saveas(gcf, sprintf('./transfer_variational/q1_%f_endtime_%f.jpg', T, T_end));
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
            saveas(gcf, sprintf('./transfer_variational/q2_%f_endtime_%f.jpg', T, T_end))
            %%%VIDEO%%%
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
        save('transfer_variational/loss_list_vi.mat', 'loss_list');
        save('transfer_variational/T_list_VI.mat', 'T_list');
    else % Solve TO
        if ~exist('./inputs_variational/', 'dir')
           mkdir('./inputs_variational/')
        end
        for j = 1:length(T_list)
            tic
            options = optimoptions('fminunc', 'MaxFunctionEvaluations', 100000, 'display', 'off');
            T_ = T_list(j);
            U = zeros(T_end/T_list(j), 1);
            q_obj = zeros(2*T_end/T_,1); % Target
            U = fminunc(@(U)ObjFunc2(U, q_0, q_m,q_obj, T_, T_end), U,options);
            qf = forward2(q_0,q_m, U, T_, T_end);  
            loss = norm(qf(2*T_end/T_-1:2*T_end/T_) - q_obj(2*T_end/T_-1:2*T_end/T_));
            save(sprintf('./inputs_variational/U_%f.mat', T_), 'U');
            T_simul = T_;
            fprintf("TO done at h : %f\n", T_);
            fprintf("Loss : %f", loss);
            toc
        end
    end
end

%% Functions

function q = forward2(q_0,q_m, U, T, T_end)
% Forward simulation
    for k = 1:int64(T_end/T)
        if k == 1
            q(2*k-1:2*k, 1) = solve_dynamics2(q_0, q_m, [0;0], [U(k);0], T);
        elseif k == 2
            q(2*k-1:2*k, 1) = solve_dynamics2(q(2*k-3:2*k-2), q_0, [U(k-1);0], [U(k); 0], T);
        else
            q(2*k-1:2*k, 1) = solve_dynamics2(q(2*k-3:2*k-2), q(2*k-5:2*k-4), [U(k-1);0], [U(k); 0], T);
        end
    end
end

function f = ObjFunc2(U, q_0, q_m, q_obj, T, T_end) % objective function without gradient
    % Calculate objective
    % Forward simulation
    q = forward2(q_0, q_m,U, T, T_end);
    f = norm(q(2*T_end/T-1:2*T_end/T) - q_obj(2*T_end/T-1:2*T_end/T));
end

function q = solve_dynamics2(q_, q_0, u_plus, u_minus, T_)
    q = fsolve(@(x)varint2(x, q_, q_0, u_plus, u_minus, T_), q_, optimoptions('fsolve', 'display', 'off'));
  
end
function f = varint2(x, q_, q_0, u_plus, u_minus, T_) 
q1x = x(1);
q1y = x(2);
q0x = q_(1);
q0y = q_(2);
dm = [(((10*cos(q0y/2 + q1y/2) + 688/45)*(conj(q0x) - conj(q1x)))/(2*conj(T_)) + ((5*cos(q0y/2 + q1y/2) + 463/90)*(conj(q0y) - conj(q1y)))/(2*conj(T_)))/T_ - 49*sin(q0x/2 + q1x/2) - (49*sin(q0x/2 + q1x/2 + q0y/2 + q1y/2))/2 + ((10*cos(q0y/2 + q1y/2) + 688/45)*(q0x - q1x))/(2*T_*conj(T_)) + ((5*cos(q0y/2 + q1y/2) + 463/90)*(q0y - q1y))/(2*T_*conj(T_)), ((463*(conj(q0y) - conj(q1y)))/(180*conj(T_)) + ((5*cos(q0y/2 + q1y/2) + 463/90)*(conj(q0x) - conj(q1x)))/(2*conj(T_)))/T_ - (49*sin(q0x/2 + q1x/2 + q0y/2 + q1y/2))/2 - ((q0x - q1x)*((5*sin(q0y/2 + q1y/2)*(conj(q0x) - conj(q1x)))/(2*conj(T_)) - (5*cos(q0y/2 + q1y/2) + 463/90)/(2*conj(T_)) + (5*sin(q0y/2 + q1y/2)*(conj(q0y) - conj(q1y)))/(4*conj(T_))))/T_ + ((463/(180*conj(T_)) - (5*sin(q0y/2 + q1y/2)*(conj(q0x) - conj(q1x)))/(4*conj(T_)))*(q0y - q1y))/T_]; 
q1x = q_(1);
q1y = q_(2);
q0x = q_0(1);
q0y = q_0(2);
dp =[- (49*sin(q0x/2 + q1x/2 + q0y/2 + q1y/2))/2 - 49*sin(q0x/2 + q1x/2) - (((10*cos(q0y/2 + q1y/2) + 688/45)*(conj(q0x) - conj(q1x)))/(2*conj(T_)) + ((5*cos(q0y/2 + q1y/2) + 463/90)*(conj(q0y) - conj(q1y)))/(2*conj(T_)))/T_ - ((10*cos(q0y/2 + q1y/2) + 688/45)*(q0x - q1x))/(2*T_*conj(T_)) - ((5*cos(q0y/2 + q1y/2) + 463/90)*(q0y - q1y))/(2*T_*conj(T_)), - (49*sin(q0x/2 + q1x/2 + q0y/2 + q1y/2))/2 - ((463*(conj(q0y) - conj(q1y)))/(180*conj(T_)) + ((5*cos(q0y/2 + q1y/2) + 463/90)*(conj(q0x) - conj(q1x)))/(2*conj(T_)))/T_ - ((q0x - q1x)*((5*cos(q0y/2 + q1y/2) + 463/90)/(2*conj(T_)) + (5*sin(q0y/2 + q1y/2)*(conj(q0x) - conj(q1x)))/(2*conj(T_)) + (5*sin(q0y/2 + q1y/2)*(conj(q0y) - conj(q1y)))/(4*conj(T_))))/T_ - ((463/(180*conj(T_)) + (5*sin(q0y/2 + q1y/2)*(conj(q0x) - conj(q1x)))/(4*conj(T_)))*(q0y - q1y))/T_];
 
f = double(T_*(dp + dm) - (u_plus+u_minus)');
end

function m = Mass(q,r, m, L, J) % mass matrix
    m = [m(1)*r(1)^2*L(1)^2+m(2)*r(2)^2*L(2)^2+m(2)*L(1)^2+2*m(2)*r(2)*L(1)*L(2)*cos(q(2)), m(2)*(r(2)^2*L(2)^2+r(2)*L(1)*L(2)*cos(q(2)));
        m(2)*(r(2)^2*L(2)^2+r(2)*L(1)*L(2)*cos(q(2))), m(2)*r(2)^2*L(2)^2] + [J(1)+J(2), J(2); J(2), J(2)];
end

