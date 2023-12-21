% -----------------------------------------------------------------
%  number_simulatios_teste.m
% -----------------------------------------------------------------
%  programmer: Diego Matos Silva Lopes
%              diego.matos@uerj.br
%
%  last update: Nov 21, 2022
% -----------------------------------------------------------------
%  
% =================================================================

clear all  
clc
figpath = '../figures/';
addpath('./utils');
addpath('./equations');

% -----------------------------------------------------------------
disp('=============================================================================');
disp(' Duffing routine that calculate the mean of RMSE varying the number of tests ');
disp('=============================================================================');

% ===============================================================
% Rossler Equation
% x_1' = x_2
% x_2' = delta*x_2 + alpha*x_1 + beta*(x_1)^3 + gamma*cos(x_3)
% x_3' = omega
% ===============================================================

%% Parameters

% Fixing the seed of randomly generated numbers
rng(16081995)

% Parameters of the Duffing dynamics

delta =  0.1;
alpha = -1.0;
beta  =  1.0;
gamma =  2.0;
omega =  1.0;

param = [delta,alpha,beta,gamma,omega];

% Initial Condition to training the SINDy
IC_t = [1, -1, 0];

% Initial Condition to test the 
IC = [2,    0,   0;
      0,   -2,   0;
     -2,    0,   0];

% Parameters SINDy
lambda = 0.085;        % threshold parameter
polyorder = 5;        % Polynomial order in the library
usesine = 1;          % Use trigonometric terms
n = 3;                % Dimension

% Parameters of the data set
t_ini = 0.0;
dt_int = 0.01;
t_fin = 30.0;
noise = 0.001;
dt_jump = 10;

% Parameters of the MSE 
t_ini_mse = 0.0;
t_fin_mse = 100.0;
dt_int_mse = 0.1;
tspan_mse = [t_ini_mse:dt_int_mse:t_fin_mse];
repeat = [10, 25, 50, 100, 250, 500];
%repeat = [10, 25, 50];

%% Generate Data

% Integrate
tic;
for l=1:size(IC,1)

    tspan=(t_ini:dt_int:t_fin);
    total_data_size = length(tspan);
    options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
    [t,x]=ode45(@(t,x)duffing(t,x,param),tspan,IC_t,options);
    
    % Extract sparser subset
    cont = 1;
    for i=1:dt_jump:(total_data_size)
        data_t(cont,1) = t(i,1);
        data_x(cont,:) = x(i,:);
        cont = cont + 1;
    end
    
    for k=1:length(repeat)
        for j=1:repeat(k)
            % Clear the matrix data
            data_dx = [];
            xA = []; tA = [];
            xB = []; tB = [];
    
            %% Compute Derivative
            for i=1:length(data_x)
                data_dx(i,:) = duffing(0,data_x(i,:),param);
            end
            data_dx = data_dx + noise*randn(size(data_dx));
 
            %% Pool Data  (i.e., build library of nonlinear time series)

            Theta = poolData(data_x,n,polyorder,usesine);

            %% Compute Sparse Regression: sequential least squares
         
            %disp('========== Result of the SINDy with STLS sparse promoter ===========');
            Xi = sparsifyDynamics(Theta,data_dx,lambda,n);
            %disp(Xi)

            %% Calculate the MSE
    
            % Integrate True and Identified Systems for the correlation
            [tA,xA] = ode45(@(t,x)duffing(t,x,param),tspan_mse,IC(l,:),options);   % true model
            [tB,xB] = ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan_mse,IC(l,:),options);  % approximate
    
            rmse_error(j,k) = sqrt(immse(xA,xB));
            
        end
    end
    disp(rmse_error)
    for i = 1:length(repeat)
        for p = 1:repeat(i)
            rmse_temp(p) = rmse_error(p,i);
        end
        mean_rmse_error(l,i) = mean(rmse_temp);
    end
    disp(mean_rmse_error)
end  

toc

figure
stem(transpose(repeat), transpose([mean_rmse_error(1,:); mean_rmse_error(2,:); mean_rmse_error(3,:)]),'LineStyle','--')
%plot(repeat,mean_rmse_error(1,:),'g-o',repeat,mean_rmse_error(2,:),'b-o',repeat,mean_rmse_error(3,:),'r-o')
legend({'$IC_1$','$IC_2$','$IC_3$'},'Interpreter','latex','fontsize',18)
xlabel('Number of time series measurements','FontSize',13)
ylabel('RMSE','FontSize',13)
ylim([0 0.08])
xlim([0 510])
set(gcf,'Position',[400 175 400 175])
print(gcf,'Sim_complete_100_lim2','-dpdf','-r300','-bestfit')
system(['pdfcrop',' ','N_sim_complete_100_lim','.pdf',' ','N_sim_complete_100_lim','.pdf'])