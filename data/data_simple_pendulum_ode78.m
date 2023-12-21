% Generate and save pendulum data
% Diego Matos S. L.
% last update 08/06/20
% UERJ

clear all, close all, clc

%% Generate Data

% Parameters
rng(16081995)    % fixing the seed of randomly generated numbers
n = 2;           % dimension
ang = 1.0;       % initial position in radians
v_ang = 1.0;     % initial velocity in radians/s
g = 9.81;        % gravity
l = 9.81;        % length line
ti = 0;          % initial time
tf = 5;         % final time
int_dt = .01;    % time step fot the integration
dt_jump = 10;    % jump in the time step for the sparse dataset
eps = 0.0;       % noise intensity

param = [g; l];  % parameters equation
x0 = [ang, v_ang];     % initial condition

%% Integrate

tspan=[ti:int_dt:tf];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,x]=BV78(@(t,x)simple_pendulum(t,x,param),tspan,x0,1e-12,1e-12*ones(1,n),false);

% Sparse Dataset
cont = 1;
for i=1:dt_jump:N
    data_t(cont,1) = t(i,1);
    data_theta(cont,:) = x(i,:);
    cont = cont+1;
end

%% Compute Derivative

for i=1:length(data_theta)
    data_dtheta(i,:) = simple_pendulum(data_t(i,1),data_theta(i,:),param);
end
data_dtheta = data_dtheta + eps*randn(size(data_dtheta));

%% Saving the data

save dtheta_simple_pendulum_ode78_5.dat data_dtheta -ascii
save theta_simple_pendulum_ode78_5.dat data_theta -ascii
save t_simple_pendulum_ode78_5.dat data_t -ascii