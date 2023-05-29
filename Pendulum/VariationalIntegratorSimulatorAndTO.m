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

T_list = [0.2 0.1 0.05 0.02 0.01 0.008 0.005]; % step sizes to use for TO
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


% For simulation only
T_simul = 0.0001; % Time step size of main simulation
simul = true; % false means solving TO to get optimal inputs under step sizes of T_list

%% Initialize
% initial conditions, q(-1) and q(0)
q_0 = [pi/2;0]; % initial condition q(0)
q_m = [pi/2;0]; % initial condition q(-1) -> v(0) = 0

%% TO

if simul == true % uses inputs calculated from TO and transfer it to model under step size of T_simulation
    if ~exist('./transfer_variational/', 'dir')
       mkdir('./transfer_variational/')
    end
    for j = 1:length(T_list)
        T_ = T_list(j);
        load(sprintf('./inputs_variational/U_%f.mat',T_), 'U');
        for i = 1:T_end/T_
            U_simul((T_ / T_simul)*(i-1)+1:(T_ / T_simul)*i,1) = U(i) * (T_simul / T_);
            %%% input is scaled by (T_simul / T_)
        end
        qf = forward2(q_0, q_m, U_simul, T_simul, T_end);
        q_obj = zeros(2*T_end/T_simul,1); % Target
        loss = norm(qf(2*T_end/T_simul-1:2*T_end/T_simul) - q_obj(2*T_end/T_simul-1:2*T_end/T_simul));
        loss_list = [loss_list, loss];
        save(sprintf('./transfer_variational/q_%f_simultime_%f.mat', T_, T_simul), 'qf');
        save(sprintf('./transfer_variational/loss_%f_simultime_%f.mat', T_, T_simul), 'loss');
        fprintf("Simulation complete at h : %f\n", T_);
        fprintf("Loss : %f\n", loss)
    end
    save('./transfer_variational/loss_list_VI', 'loss_list');
    save('./transfer_variational/T_list_VI', 'T_list');
else % TO MODE
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


%% Plot
% T_ = T_simul;
% t = T_:T_:T_end;
% plot(t, qf(1:2:2*T_end/T_-1));
% hold on;
% plot(t, q_obj(1:2:2*T_end/T_-1));
% xlabel('time, [s]');
% ylabel('angle, [rad]');
% title('q_1');
% legend('real', 'target');
% figure(2)
% plot(t, qf(2:2:2*T_end/T_));
% hold on;
% plot(t, q_obj(2:2:2*T_end/T_));
% xlabel('time, [s]');
% ylabel('angle, [rad]');
% title('q_2');
% legend('real', 'target');
% 
% v = VideoWriter('TO_variational.avi');
% open(v);
% for k = T_end/T_/50:T_end/T_/50:T_end/T_
%     figure(3)
%     q1 = qf(2*k-1);
%     q2 = qf(2*k);
%     x = zeros(2, 3);
%     x(:, 2) = L_(1)*[sin(q1); -cos(q1)];
%     x(:, 3) = x(:, 2) + L_(2)*[sin(q1+q2); -cos(q1+q2)];
% 
%     clf;
%     plot(x(1, :), x(2, :), 'o-');
%     axis equal;
%     xlim([-1, 1]);
%     ylim([-2.5,2.5]);
%     drawnow;
%    frame = getframe(gcf);
%    writeVideo(v,frame);
% end
% close(v);

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

