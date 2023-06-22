clear all; close all; clc;
%% Simulation Settings
N = 2; % number of bodies
T_end = 10; % simulation end time

T = 0.001;

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

q_0 = [pi/2;0]; % initial condition q(0)
q_m = [pi/2;0]; % initial condition q(-1) -> v(0) = 0
M = zeros(2,2*T_end/T); % Mass matrix
C = zeros(2,2*T_end/T); % C&C matrix 
G = zeros(2,T_end/T); % Gravity vector
U = zeros(T_end/T,1); % Torque input vector. [u1]

% initial conditions, q(-1) and q(0)
v_0 = [0;0]; % initial condition q(-1) -> v(0) = 0


% For simulation only
[q_euler, ~, ~, ~, E_euler] = forward_euler(q_0,q_m, r, U, m, L, g, J, T, T_end);
[q_pmi, v_pmi, ~, ~, ~, E_pmi] = forward_pmi(q_0, v_0, r, U, m, L, g, J, T, T_end);
[q_vi, E_vi] = forward_vi(q_0,q_m, U, r, m, L, J, g, T, T_end);
figure(1);
plot(T:T:T_end, E_euler);
hold on;
plot(T:T:T_end, E_pmi);
hold on;
plot(T:T:T_end, E_vi);
xlabel('time, [s]');
ylabel('Energy, [j]');
title('Energy property');
legend('Euler', 'PMI', 'VI');

for k = 200 :200 :T_end/T
    figure(2)
    q1 = q_vi(2*k-1);
    q2 = q_vi(2*k);
    x = zeros(2, 3);
    x(:, 1) = [0;0];
    x(:, 2) = [L(1)*sin(q1) ; -L(1)*cos(q1)];
    x(:, 3) = x(:, 2) + L(2)*[sin(q1+q2); -cos(q1+q2)];
    
    clf;
    plot(x(1, :), x(2, :), 'o-');
    axis equal;
    xlim([-3, 3]);
    ylim([-3,3]);
    drawnow;
    frame = getframe(gcf);
end

%%%%%%%%%%%FUNCTIONS%%%%%%%%%%%%%
function [q, v, M, C, G, E] = forward_pmi(q_0,v_0, r, U, m, L, g, J, T, T_end)
    % Forward simulation
    q = zeros(2*int64(T_end/T),1);
    v = zeros(2*int64(T_end/T),1);
    for k = 1:int64(T_end/T)
        if (k == 1)
            [q(2*k-1:2*k), v(2*k-1:2*k)] = solve_dynamics_pmi(q_0, v_0,r, [U(k);0], m, L, g, J, T);
        else
            [q(2*k-1:2*k), v(2*k-1:2*k)] = solve_dynamics_pmi(q(2*k-3:2*k-2), v(2*k-3:2*k-2),r, [U(k);0], m, L, g, J, T);
        end
    % Calculate dynamics matrices
        M(:, 2*k-1:2*k) = Mass(q(2*k-1:2*k),r, m, L, J);
        C(:, 2*k-1:2*k) = Cori(q(2*k-1:2*k), v(2*k-1:2*k) ,r,  m, L);
        G(:, k) = Grav(q(2*k-1:2*k),r, m, L, g);
        E(k) = 0.5*v(2*k-1:2*k)'*M(:, 2*k-1:2*k)*v(2*k-1:2*k) - (L(1)*r(1)*m(1)*g*cos(q(2*k-1))+(L(1)*r(1)*cos(q(2*k-1))+L(2)*r(2)*cos(q(2*k-1)+q(2*k)))*m(2)*g);
    end
end

function [q,v] = solve_dynamics_pmi(q_0_, v_0_,r_, u_, m_, L_, g_, J_, T_) % Forward simulation works well
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


function [q, M, C, G, E] = forward_euler(q_0,q_m, r, U, m, L, g, J, T, T_end)
            % Forward simulation
    q = zeros(2*int64(T_end/T),1);
    for k = 1:int64(T_end/T)
        if k == 1
            q(2*k-1:2*k) = solve_dynamics_euler(q_0, q_m,r, [U(k);0], m, L, g, J, T);
            v = (q(2*k-1:2*k) - q_0)/T;
        elseif k == 2
            q(2*k-1:2*k) = solve_dynamics_euler(q(2*k-3:2*k-2), q_0, r, [U(k);0], m, L, g, J, T);
            v = (q(2*k-1:2*k) - q(2*k-3:2*k-2))/T;
        else
            q(2*k-1:2*k) = solve_dynamics_euler(q(2*k-3:2*k-2), q(2*k-5:2*k-4),r, [U(k);0], m, L, g,J, T);            
            v = (q(2*k-1:2*k) - q(2*k-3:2*k-2))/T;
        end
    % Calculate dynamics matrices
        if k == 1
            C(:, 2*k-1:2*k) = Cori(q(2*k-1:2*k), (q(2*k-1:2*k)-q_0)/T,r,  m, L);
        else
            C(:, 2*k-1:2*k) = Cori(q(2*k-1:2*k), (q(2*k-1:2*k) - q(2*k-3:2*k-2))/T, r, m, L);
        end 
        G(:, k) = Grav(q(2*k-1:2*k),r, m, L, g);
        M(:, 2*k-1:2*k) = Mass(q(2*k-1:2*k),r, m, L, J);
        E(k) = 0.5*v'*M(:, 2*k-1:2*k)*v - (L(1)*r(1)*m(1)*g*cos(q(2*k-1))+(L(1)*r(1)*cos(q(2*k-1))+L(2)*r(2)*cos(q(2*k-1)+q(2*k)))*m(2)*g);

    end
end

function solve = solve_dynamics_euler(q_0_, q_m_,r_, u_, m_, L_, g_, J_, T_) % Forward simulation works well
    m = Mass(q_0_,r_, m_, L_, J_);
    c = Cori(q_0_, (q_0_-q_m_)/T_,r_, m_, L_);
    dg = dg_dq(q_0_, r_,m_, L_, g_);
    solve = q_0_ + (m/T_+c)\(u_+m*(q_0_-q_m_)/T_^2 - Grav(q_0_, r_, m_, L_, g_))*T_;
end

function [q, E] = forward_vi(q_0,q_m, U,r, m, L, J, g, T, T_end)
% Forward simulation
    for k = 1:int64(T_end/T)
        if k == 1
            q(2*k-1:2*k, 1) = solve_dynamics_vi(q_0, q_m, [0;0], [U(k);0], T);
            qrep = (q(2*k-1:2*k)+q_0)/2;
            M = Mass(qrep,r, m, L, J);
            v = (q(2*k-1:2*k) - q_0)/T;
            V = - (L(1)*r(1)*m(1)*g*cos(qrep(1))+(L(1)*r(1)*cos(qrep(1))+L(2)*r(2)*cos(qrep(1)+qrep(2)))*m(2)*g);

        elseif k == 2
            q(2*k-1:2*k, 1) = solve_dynamics_vi(q(2*k-3:2*k-2), q_0, [U(k-1);0], [U(k); 0], T);
            qrep = (q(2*k-1:2*k)+q(2*k-3:2*k-2))/2;
            M = Mass(qrep,r, m, L, J);
            v = (q(2*k-1:2*k) - q(2*k-3:2*k-2))/T;
            V = - (L(1)*r(1)*m(1)*g*cos(qrep(1))+(L(1)*r(1)*cos(qrep(1))+L(2)*r(2)*cos(qrep(1)+qrep(2)))*m(2)*g);
        else
            q(2*k-1:2*k, 1) = solve_dynamics_vi(q(2*k-3:2*k-2), q(2*k-5:2*k-4), [U(k-1);0], [U(k); 0], T);
            qrep = (q(2*k-1:2*k)+q(2*k-3:2*k-2))/2;
            M = Mass(qrep,r, m, L, J);
            v = (q(2*k-1:2*k) - q(2*k-3:2*k-2))/T;
            V = - (L(1)*r(1)*m(1)*g*cos(qrep(1))+(L(1)*r(1)*cos(qrep(1))+L(2)*r(2)*cos(qrep(1)+qrep(2)))*m(2)*g);
        end
        E(k) = 0.5*v'*M *v +V;
    end
end

function q = solve_dynamics_vi(q_, q_0, u_plus, u_minus, T_)
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
