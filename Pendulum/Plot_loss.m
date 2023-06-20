clear all; close all; clc;
load('./transfer_euler/T_list_Euler.mat', 'T_list');
load('./transfer_euler/loss_list_Euler.mat', 'loss_list');
T_list_E = T_list;
loss_list_E = loss_list;
load('./transfer_variational/T_list_VI.mat', 'T_list');
load('./transfer_variational/loss_list_VI.mat', 'loss_list');
T_list_V = T_list;
loss_list_V = loss_list;
load('./transfer_pmi/T_list_PMI.mat', 'T_list');
load('./transfer_pmi/loss_list_PMI.mat', 'loss_list');
T_list_P = T_list;
loss_list_P = loss_list;

semilogx(T_list_V,loss_list_V, 'linewidth', 3);
hold on;
semilogx(T_list_E, loss_list_E, 'linewidth', 3);
hold on;
semilogx(T_list_P, loss_list_P, 'linewidth', 3);
legend('variational', 'Euler', 'PMI');
xlabel('time step size, [s]');
ylabel('loss');
hold on;
title('Input transfer in double pendulum problem');
ylim([0, 0.01])
