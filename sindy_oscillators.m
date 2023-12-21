%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Diego Matos S. L.                   %
% Last update 04/01/21                        %
% Master's degree student at UERJ             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
addpath('./equations');
addpath('./utils');

%% Harmonic oscillator
% x" = -(c/m)*x' - (k/m)*x + (f/m)*cos(w*t)

%m = 1.0;    % oscillator mass
%c = 0.1;    % damping term
%k = 1.0;    % spring stiffness
%f = 1.0;    % harmonic force intensity
%w = 2.0;    % frequency in the cossine harmonic force

%param = [c/m; k/m; f/m; w];
%x0 = [0.0; 2.0; -2.0]; % [initial frequency*time; displacement; velocity]

%equation = @(t,x)harmonic_oscillator(t,x,param);

%% Duffing oscillator
% x" = -d*x' - b*x - a*x^3 + g*cos(w*t)

% Without chaos
c = 0.1;          % damping
a = -1.0;         % linear stiffness
b = 1.0;          % non-linear stiffness
g = 1.0;          % harmonic force intensity
w = 2.0;          % frequency in the cossine harmonic force

param = [c, a, b, g, w];
x0 = [2.0; -2.0; 0.0]; % [displacement; velocity, initial frequency*time]

% Chaos
%c = 0.1;          % damping
%a = -1.0;         % linear stiffness
%b = 0.25;         % non-linear stiffness
%g = 2.5;          % harmonic force intensity
%w = 2.0;          % frequency in the cossine harmonic force

%param = [c, a, b, g, w];
%x0 = [1.0; -1.0; 0.0]; % [displacement; velocity, initial frequency*time]

equation = @(t,x)duffing(t,x,param);

%% Van der Pol oscillator
% x" = u(1-x)*x' - a*x + f*cos(w*t)

%u = 4.0;           % non-linear term
%a = 1.0;           % generic term
%f = 0.0;           % harmonic forve intensity
%w = 0.0;           % frequency in the cossine harmonic force

%param= [u, a, f, w];
%x0 = [0.5; -0.2; 0.0]; % [displacement; velocity; initial frequency*time]

%equation = @(t,x)van_der_pol(t,x,param);

%% Parameters SINDy

lambda = 0.05;        % threshold parameter
polyorder = 5;        % polynomial order in the library
usesine = 1;          % Use trigonometric terms
n = 3;                % dimension of the identified equation

%% Data generated

% Parameters
rng(16081995)         % fixing the seed of randomly generated numbers
ti = 0;               % initial time of the data
tf = 15;              % final time of the data
int_dt = .01;         % time step fot the integration
dt_jump = 10;         % jump in the int_dt for the sparse dataset
noise = 0.01;         % noise intensity

% Integrate the equation
tspan=[ti:int_dt:tf];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,x]=ode45(@(t,x)duffing(t,x,param),tspan,x0,options);
%[t,x]=ode45(@(t,x)van_der_pol(t,x,param),tspan,x0,options);

% Sparse Dataset
cont = 1;
for i=1:dt_jump:(N+1)
    data_t(cont,1) = t(i,1);
    data_x(cont,:) = x(i,:);
    cont = cont+1;
end

%% Compute derivative

for i=1:length(data_x)
    dx(i,:) = duffing(data_t(i,1),data_x(i,:),param);
    %dx(i,:) = van_der_pol(data_t(i,1),data_x(i,:),param);
end
dx = dx + noise*randn(size(dx));

%% Pool data (build library of nonlinear time series)

Theta = poolData(data_x,n,polyorder,usesine);

%% Compute sparse regression: STLS

Xi = sparsifyDynamics(Theta,dx,lambda,n)
%poolDataLIST({'x_1','x_2','x_3'},Xi,n,polyorder,usesine);

%% Validation

%x0 = [-3.0; 1.0; 0.0];  % Test another initial condition for the identified Ev. Law 

%% Time for plot and compare the results
ti_p = 0;
tf_p = 60;
dt_p = 0.01;
tspan = [ti_p:dt_p:tf_p];

%% Integrate true and identified systems

[t_true,x_true]=ode45(equation,tspan,x0,options);            % true model
[t_iden,x_iden]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);  % identified

%% Plots time series
if tf_p <= 0
    error('Input must be positive')
end
x_order = floor(log10(tf_p)) + 1;
y_order_dis = floor(log10(max(max(x_true(:,1),x_iden(:,1))))) + 1;
y_order_vel = floor(log10(max(max(x_true(:,2),x_iden(:,2))))) + 1;
y_max_dis = ceil(max(max(x_true(:,1),x_iden(:,1))) * 1.7);
y_min_dis = floor(min(min(x_true(:,1),x_iden(:,1))) * 1.7);
y_max_vel = ceil(max(max(x_true(:,2),x_iden(:,2))) * 1.7);
y_min_vel = floor(min(min(x_true(:,2),x_iden(:,2))) * 1.7);

% Automatic plot
% Position
plot_dynamics(t_true,x_true(:,1),t_iden,x_iden(:,1),data_t,data_x(:,1),'','Time','x_1','time_series_x1','r',x_order,y_order_dis,ti_p,tf_p,y_min_dis,y_max_dis);
% Velocity
plot_dynamics(t_true,x_true(:,2),t_iden,x_iden(:,2),data_t,data_x(:,2),'','Time','x_2','time_series_x2','b',x_order,y_order_vel,ti_p,tf_p,y_min_vel,y_max_vel);
% Phase space
%plot_dynamics(x_true(:,1),x_true(:,2),x_iden(:,1),x_iden(:,2),data_x(:,1),data_x(:,2),'','Displacement','Velocity','phase_space','c',y_order_dis,y_order_vel,y_min_dis,y_max_dis,y_min_vel,y_max_vel);

% 3d Plot
%plot_3d_ph(x_true(:,1),x_true(:,2),t_true,x_iden(:,1),x_iden(:,2),t_iden,data_x(:,1),data_x(:,2),data_t,"x_1","x_2","x_3","","test");

% Manual plot
% % Position
%plot_dynamics(t_true,x_true(:,1),t_iden,x_iden(:,1),data_t,data_x(:,1),'','Time','Displacement','disp_deslocado','r',3,1,745,755,-3,3);
% % Velocity
%plot_dynamics(t_true,x_true(:,2),t_iden,x_iden(:,2),data_t,data_x(:,2),'','Time','Velocity','vel_deslocado','b',3,1,745,755,-5,5);
% % Phase space
%plot_dynamics(x_true(:,2),x_true(:,3),x_iden(:,2),x_iden(:,3),data_x(:,2),data_x(:,3),'Phase_space','','','','c',1,1,-2.5,2.5,-4.5,4.5);