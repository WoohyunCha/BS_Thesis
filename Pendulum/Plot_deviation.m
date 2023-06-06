clear all; close all; clc;

figure(1);
T_list = [0.1 0.05 0.02 0.01 0.008 0.004 0.002];
T_end = 1;
dev_list = zeros(length(T_list), 3);
for j = 1:length(T_list)
    load(sprintf('./test_variational/q_%f_endtime_%f.mat', T_list(j), T_end), sprintf('q_%d', j));
end
for j = 2:length(T_list)
    q = eval(sprintf('q_%d', j));
    q1 = q(1:2:end);
    q2 = q(2:2:end);
    q_bench = eval(sprintf('q_%d', length(T_list)));
    dev_list(j, 1) = norm(q_bench(1:2*T_list(j)/T_list(end):end) - q1)^2 ...
                + norm(q_bench(2:2*T_list(j)/T_list(end):end) - q2)^2;
end
semilogx(T_list(2:end), dev_list(2:end, 1));
hold on;
for j = 1:length(T_list)
    load(sprintf('./test_pmi/q_%f_endtime_%f.mat', T_list(j), T_end), sprintf('q_%d', j));
end
for j = 1:length(T_list)
    q = eval(sprintf('q_%d', j));
    q1 = q(1:2:end);
    q2 = q(2:2:end);
    q_bench = eval(sprintf('q_%d', length(T_list)));
    dev_list(j, 2) = norm(q_bench(1:2*T_list(j)/T_list(end):end) - q1)^2 ...
                + norm(q_bench(2:2*T_list(j)/T_list(end):end) - q2)^2;
end
semilogx(T_list(1:end), dev_list(:, 2));
for j = 1:length(T_list)
    load(sprintf('./test_euler/q_%f_endtime_%f.mat', T_list(j), T_end), sprintf('q_%d', j));
end
for j = 1:length(T_list)
    q = eval(sprintf('q_%d', j));
    q1 = q(1:2:end);
    q2 = q(2:2:end);
    q_bench = eval(sprintf('q_%d', length(T_list)));
    dev_list(j, 3) = norm(q_bench(1:2*T_list(j)/T_list(end):end) - q1)^2 ...
                + norm(q_bench(2:2*T_list(j)/T_list(end):end) - q2)^2;
end
semilogx(T_list(1:end), dev_list(:, 3));
legend('variational integrator', 'PMI', 'Euler');
title('Deviation from benchmark trajectory');
saveas(gcf, 'plot_deviation.jpg');