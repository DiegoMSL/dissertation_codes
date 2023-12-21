% -----------------------------------------------------------------
%  duffing__intensity_noise.m
% -----------------------------------------------------------------
%  programmer: Diego Matos Silva Lopes
%              diego.matos@uerj.br
%
%  last update: jun 06, 2022
% -----------------------------------------------------------------
%  
% =================================================================

clear all  
clc
figpath = '../figures/';
addpath('./utils');
addpath('./equations')

% -----------------------------------------------------------------
disp('=================================================================');
disp(' Duffing routine that calculate the mean of MSE varying the noise');
disp('=================================================================');

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

delta = 0.1;
alpha = -1.0;
beta = 1.0;
gamma = 2.0;
omega = 2.0;
param = [delta,alpha,beta,gamma,omega];

% Initial Condition
IC_t = [1, -1, 0];

IC = [1,    0,   0;
      1.5, -0.5, 0;
      2,   -1,   0];
      %-1,   0,   0;
      %-2,   1,   0];

% Parameters SINDy
lambda = 0.085;        % threshold parameter
polyorder = 4;        % Polynomial order in the library
usesine = 1;          % Use trigonometric terms
n = 3;                % Dimension

% Parameters od the data set
t_ini = 0.0;
dt_int = 0.001;
t_fin = 30.0;
noise = [0, 0.005, 0.010, 0.015, 0.025, 0.035, 0.050];
dt_jump = 50;

% Parameters of the MSE 
t_ini_mse = 0.0;
t_fin_mse = 100.0;
dt_int_mse = 0.1;
tspan_mse = (t_ini_mse:dt_int_mse:t_fin_mse);
repeat = 100;

% Parameters of the correlation test 
t_ini_corr = 0.0;
t_fin_corr = 100.0;
dt_int_corr = 0.1;
tspan_corr = [t_ini_corr:dt_int_corr:t_fin_corr];

%% Generate Data

% Integrate

for l=1:size(IC,1)
    
    tspan=(t_ini:dt_int:t_fin);
    total_data_size = length(tspan);
    options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
    [t,x]=ode45(@(t,x)duffing(t,x,param),tspan,IC_t,options);
    
    for k=1:repeat    
        for j=1:length(noise)
            % Clear the matrix data
            data_t = [];
            data_x = [];
            data_dx = [];
            xA = []; tA = [];
            xB = []; tB = [];
    
            % Extract sparser subset
            cont = 1;
            for i=1:dt_jump:(total_data_size)
                data_t(cont,1) = t(i,1);
                data_x(cont,:) = x(i,:);
                cont = cont+1;
            end
        
            %% Compute Derivative
            for i=1:length(data_x)
                data_dx(i,:) = duffing(0,data_x(i,:),param);
            end
            
            data_dx = data_dx + noise(j)*randn(size(data_dx));
            
            text_noise = ['intensity of noise used = ',num2str(noise(j))];
            disp(text_noise)
            
            %% Pool Data  (i.e., build library of nonlinear time series)

            Theta = poolData(data_x,n,polyorder,usesine);

            %% Compute Sparse Regression: sequential least squares
         
            disp('========== Result of the SINDy with STLS sparse promoter ===========');
            Xi = sparsifyDynamics(Theta,data_dx,lambda,n);
            disp(Xi)

            %% Calculate the MSE
    
            % Integrate True and Identified Systems for the correlation
            [tA,xA] = ode45(@(t,x)duffing(t,x,param),tspan_mse,IC(l,:),options);   % true model
            [tB,xB] = ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan_mse,IC(l,:),options);  % approximate
    
            rmse_error(k,j) = sqrt(immse(xA,xB));
            
            % Integrate True and Identified Systems for the correlation
            [tA_corr,xA_corr] = ode45(@(t,x)duffing(t,x,param),tspan_corr,IC(l,:),options);   % true model
            [tB_corr,xB_corr] = ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan_corr,IC(l,:),options);  % approximate
            
            R1 = corrcoef(xA_corr(:,1),xB_corr(:,1));
            R2 = corrcoef(xA_corr(:,2),xB_corr(:,2));
            R3 = corrcoef(xA_corr(:,3),xB_corr(:,3));
            Rt = corrcoef(xA_corr,xB_corr);
        
            R1_plot(k,j) = R1(1,2);
            R2_plot(k,j) = R2(1,2);
            R3_plot(k,j) = R3(1,2);
            Rt_plot(k,j) = Rt(1,2);

        end
    end
    
    mean_rmse_error(l,:) = mean(rmse_error);
    mean_R1_plot(l,:) = mean(R1_plot);
    mean_R2_plot(l,:) = mean(R2_plot);
    mean_R3_plot(l,:) = mean(R3_plot);
    mean_Rt_plot(l,:) = mean(Rt_plot);

end

%csvwrite('table/rmse_noise_poly_4.csv',mean_rmse_error)
%csvwrite('table/R1_noise_poly_4.csv',mean_R1_plot)
%csvwrite('table/R2_noise_poly_4.csv',mean_R2_plot)
%csvwrite('table/R3_noise_poly_4.csv',mean_R3_plot)
%csvwrite('table/Rt_noise_poly_4.csv',mean_Rt_plot)

figure
%plot(noise,mean_rmse_error(1,:),'g-o',noise,mean_rmse_error(2,:),'b-o',noise,mean_rmse_error(3,:),'c-o',noise,mean_rmse_error(4,:),'m-o',noise,mean_rmse_error(5,:),'r-o')
stem(transpose(noise), transpose([mean_rmse_error]),'LineStyle','--')
legend({'$IC_1$','$IC_2$','$IC_3$'},'Interpreter','latex','fontsize',18,'southeast')
xlabel('Noise intensity','FontSize',13)
ylabel('RMSE','FontSize',13)
set(gcf,'Position',[900 400 900 400])
print(gcf,'RMSE_mean_noise_poly_4','-dpdf','-r300','-bestfit');
system(['pdfcrop',' ','RMSE_mean_noise_poly_4','.pdf',' ','RMSE_mean_noise_poly_4','.pdf'])

figure
%plot(noise,mean_R1_plot(1,:),'g-o',noise,mean_R1_plot(2,:),'b-o',noise,mean_R1_plot(3,:),'c-o',noise,mean_R1_plot(4,:),'m-o',noise,mean_R1_plot(5,:),'r-o')
stem(transpose(noise), transpose([mean_R1_plot]),'LineStyle','--')
legend({'$IC_1$','$IC_2$','$IC_3$'},'Interpreter','latex','fontsize',18,'Location','southeast')
xlabel('Noise intensity','FontSize',13)
ylabel('Correlation $\dot{x}_1$','Interpreter','latex','FontSize',13)
set(gcf,'Position',[900 400 900 400])
print(gcf,'corr_x_mean_noise_poly_4','-dpdf','-r300','-bestfit');
system(['pdfcrop',' ','corr_x_mean_noise_poly_4','.pdf',' ','corr_x_mean_noise_poly_4','.pdf'])

figure
%plot(noise,mean_R2_plot(1,:),'g-o',noise,mean_R2_plot(2,:),'b-o',noise,mean_R2_plot(3,:),'c-o',noise,mean_R2_plot(4,:),'m-o',noise,mean_R2_plot(5,:),'r-o')
stem(transpose(noise), transpose([mean_R2_plot]),'LineStyle','--')
%legend({'$IC_1$','$IC_2$','$IC_3$'},'Interpreter','latex','fontsize',18,'Location','southeast')
xlabel('Noise intensity','FontSize',13)
ylabel('Correlation $\dot{x}_2$','Interpreter','latex','FontSize',13)
set(gcf,'Position',[900 400 900 400])
print(gcf,'corr_y_mean_noise_poly_4','-dpdf','-r300','-bestfit');
system(['pdfcrop',' ','corr_y_mean_noise_poly_4','.pdf',' ','corr_y_mean_noise_poly_4','.pdf'])

figure
%plot(noise,mean_R3_plot(1,:),'g-o',noise,mean_R3_plot(2,:),'b-o',noise,mean_R3_plot(3,:),'c-o',noise,mean_R3_plot(4,:),'m-o',noise,mean_R3_plot(5,:),'r-o')
stem(transpose(noise), transpose([mean_R3_plot]),'LineStyle','--')
legend({'$IC_1$','$IC_2$','$IC_3$'},'Interpreter','latex','fontsize',18,'Location','southeast')
xlabel('Noise intensity','FontSize',13)
ylabel('Correlation $\dot{x}_3$','Interpreter','latex','FontSize',13)
set(gcf,'Position',[900 400 900 400])
print(gcf,'corr_z_mean_noise_poly_4','-dpdf','-r300','-bestfit')
system(['pdfcrop',' ','corr_z_mean_noise_poly_4','.pdf',' ','corr_z_mean_noise_poly_4','.pdf'])

figure
%plot(noise,mean_Rt_plot(1,:),'g-o',noise,mean_Rt_plot(2,:),'b-o',noise,mean_Rt_plot(3,:),'c-o',noise,mean_Rt_plot(4,:),'m-o',noise,mean_Rt_plot(5,:),'r-o')
stem(transpose(noise), transpose([mean_Rt_plot]),'LineStyle','--')
legend({'$IC_1$','$IC_2$','$IC_3$'},'Interpreter','latex','fontsize',18,'Location','southeast')
xlabel('Noise intensity','FontSize',13)
ylabel('Phase space correlation','FontSize',13)
set(gcf,'Position',[900 400 900 400])
print(gcf,'corr_ph_mean_noise_poly_4','-dpdf','-r300','-bestfit')
system(['pdfcrop',' ','corr_ph_mean_noise_poly_4','.pdf',' ','corr_ph_mean_noise_poly_4','.pdf'])
