% Generate and save pendulum data
% Diego Matos S. L.
% last update 08/06/20
% UERJ

clear all, close all, clc
addpath('./utils');
addpath('./data');
addpath('./equations');

%% Load the data

[filename directory_name] = uigetfile('data/t_simple_pendulum_ode45_15.dat');
fullname = fullfile(directory_name, filename);
t = load(fullname);

[filename directory_name] = uigetfile('data/theta_simple_pendulum_ode45_15.dat');
fullname = fullfile(directory_name, filename);
theta = load(fullname);

[filename directory_name] = uigetfile('data/dtheta_simple_pendulum_ode45_15.dat');
fullname = fullfile(directory_name, filename);
dtheta = load(fullname);

%% parameter
n = 2;               % dimension
polyorder = 5;       % polinomial order of functions
usesine = 0;         % use trignometric terms
lambda = 0.001;
ti = 0.0;            % initial time plot
tf = 30.0;           % final time plot
interval = 0.01;      % 
tspant = [ti:interval:tf];
x0 = [theta(1,1),theta(1,2)];
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
g = 9.81;        % gravity
l = 9.81;        % length line
param = [g; l];  % parameters equation

%% pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(theta,n,polyorder,usesine);
m = size(Theta,2);

%% compute Sparse regression: sequential least squares
Xi = sparsifyDynamics(Theta,dtheta,lambda,n)

%% integrate true and identified systems
[t_true,x_true]=ode45(@(t,theta)simple_pendulum(t,theta,param),tspant,x0,options);   % true model
[t_iden,x_iden]=ode45(@(t,theta)sparseGalerkin(t,theta,Xi,polyorder,usesine),tspant,x0,options);  % approximate

%% Plots time series
if tf <= 0
    error('Input must be positive')
end
x_order = floor(log10(tf)) + 1;
y_order_dis = floor(log10(max(max(x_true(:,1),x_iden(:,1))))) + 1;
y_order_vel = floor(log10(max(max(x_true(:,2),x_iden(:,2))))) + 1;
y_max_dis = ceil(max(max(x_true(:,1),x_iden(:,1))) * 1.7);
y_min_dis = floor(min(min(x_true(:,1),x_iden(:,1))) * 1.7);
y_max_vel = ceil(max(max(x_true(:,2),x_iden(:,2))) * 1.7);
y_min_vel = floor(min(min(x_true(:,2),x_iden(:,2))) * 1.7);

% Automatic plot
% Position
plot_dynamics(t_true,x_true(:,1),t_iden,x_iden(:,1),t,theta(:,1),'','Time','Angle','time_series_angle','r',x_order,y_order_dis,ti,tf,y_min_dis,y_max_dis);
% Velocity
plot_dynamics(t_true,x_true(:,2),t_iden,x_iden(:,2),t,theta(:,2),'','Time','Angular Vel.','time_series_velocity_angle','b',x_order,y_order_vel,ti,tf,y_min_vel,y_max_vel);
% Phase space
plot_dynamics(x_true(:,1),x_true(:,2),x_iden(:,1),x_iden(:,2),theta(:,1),theta(:,2),'','Angle','Angular Velocity','phase_space','c',y_order_dis,y_order_vel,y_min_dis,y_max_dis,y_min_vel,y_max_vel);
