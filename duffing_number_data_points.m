% -----------------------------------------------------------------
%  duffing_number_data_points.m
% -----------------------------------------------------------------
%  programmer: Diego Matos Silva Lopes
%              diego.matos@uerj.br
%
%  last update: April 30, 2022
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
disp(' Duffing routine that calculate the mean of MSE varying the data ');
disp('=================================================================');

% ===============================================================
% Duffing Equation
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
beta =  1.0;
gamma = 2.0;
omega = 1.0;

% EX 2 % lambda 0.175 (não funciona)
%delta = -0.2;
%alpha = 1.0;
%beta = -1.0;
%gamma = 1.2;
%omega = 1.0;

param = [delta,alpha,beta,gamma,omega];

% Initial Condition
IC_t = [1, -1, 0];

IC = [1,    0,   0;
      1.5, -0.5, 0;
      %2,   -1,   0;
      %-1,   0,   0;
      -2,   1,   0];

% Parameters SINDy
lambda = 0.085;        % threshold parameter
polyorder = 3;        % Polynomial order in the library
usesine = 0;          % Use trigonometric terms
n = 3;                % Dimension

% Parameters of the data set
t_ini = 0.0;
dt_int = 0.001;
t_fin = 30.0;
noise = 0.005;
dt_jump = [10, 20, 30, 50, 75, 100, 150];

% Parameters of the MSE 
t_ini_mse = 0.0;
t_fin_mse = 100.0;
dt_int_mse = 0.1;
tspan_mse = [t_ini_mse:dt_int_mse:t_fin_mse];
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
        for j=1:length(dt_jump)
            % Clear the matrix data
            data_t = [];
            data_x = [];
            data_dx = [];
            xA = []; tA = [];
            xB = []; tB = [];
    
            % Extract sparser subset
            cont = 1;
            for i=1:dt_jump(j):(total_data_size)
                data_t(cont,1) = t(i,1);
                data_x(cont,:) = x(i,:);
                cont = cont+1;
            end
        
            if k == 1
                number_points(j) = length(data_x);
            end
 
            text_points = ['Number os data points used = ',num2str(length(data_x))];
            disp(text_points)
        

            %% Compute Derivative
            for i=1:length(data_x)
                data_dx(i,:) = duffing(0,data_x(i,:),param);
            end
            data_dx = data_dx + noise*randn(size(data_dx));

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

for i=1:length(number_points)
    x_plot(i) = number_points(length(number_points)-i+1);
    ic_1(i) = mean_rmse_error(1,length(number_points)-i+1);
    ic_2(i) = mean_rmse_error(2,length(number_points)-i+1);
    ic_3(i) = mean_rmse_error(3,length(number_points)-i+1);
    %ic_4(i) = mean_rmse_error(4,length(number_points)-i+1);
    %ic_5(i) = mean_rmse_error(5,length(number_points)-i+1);
        
    ic_1_R1(i) = mean_R1_plot(1,length(number_points)-i+1);
    ic_2_R1(i) = mean_R1_plot(2,length(number_points)-i+1);
    ic_3_R1(i) = mean_R1_plot(3,length(number_points)-i+1);
    %ic_4_R1(i) = mean_R1_plot(4,length(number_points)-i+1);
    %ic_5_R1(i) = mean_R1_plot(5,length(number_points)-i+1);
    
    ic_1_R2(i) = mean_R2_plot(1,length(number_points)-i+1);
    ic_2_R2(i) = mean_R2_plot(2,length(number_points)-i+1);
    ic_3_R2(i) = mean_R2_plot(3,length(number_points)-i+1);
    %ic_4_R2(i) = mean_R2_plot(4,length(number_points)-i+1);
    %ic_5_R2(i) = mean_R2_plot(5,length(number_points)-i+1);
    
    ic_1_R3(i) = mean_R3_plot(1,length(number_points)-i+1);
    ic_2_R3(i) = mean_R3_plot(2,length(number_points)-i+1);
    ic_3_R3(i) = mean_R3_plot(3,length(number_points)-i+1);
    %ic_4_R3(i) = mean_R3_plot(4,length(number_points)-i+1);
    %ic_5_R3(i) = mean_R3_plot(5,length(number_points)-i+1);
    
    ic_1_Rt(i) = mean_Rt_plot(1,length(number_points)-i+1);
    ic_2_Rt(i) = mean_Rt_plot(2,length(number_points)-i+1);
    ic_3_Rt(i) = mean_Rt_plot(3,length(number_points)-i+1);
    %ic_4_Rt(i) = mean_Rt_plot(4,length(number_points)-i+1);
    %ic_5_Rt(i) = mean_Rt_plot(5,length(number_points)-i+1); 

end    

csvwrite('table/ic_1_poly_4.csv',ic_1)
csvwrite('table/ic_2_poly_4.csv',ic_2)
csvwrite('table/ic_3_poly_4.csv',ic_3)
csvwrite('table/ic_1_R1_poly_4.csv',ic_1_R1)
csvwrite('table/ic_2_R1_poly_4.csv',ic_2_R1)
csvwrite('table/ic_3_R1_poly_4.csv',ic_3_R1)
csvwrite('table/ic_1_R2_poly_4.csv',ic_1_R2)
csvwrite('table/ic_2_R2_poly_4.csv',ic_2_R2)
csvwrite('table/ic_3_R2_poly_4.csv',ic_3_R2)
csvwrite('table/ic_1_R3_poly_4.csv',ic_1_R3)
csvwrite('table/ic_2_R3_poly_4.csv',ic_2_R3)
csvwrite('table/ic_3_R3_poly_4.csv',ic_3_R3)
csvwrite('table/ic_1_Rt_poly_4.csv',ic_1_Rt)
csvwrite('table/ic_2_Rt_poly_4.csv',ic_2_Rt)
csvwrite('table/ic_3_Rt_poly_4.csv',ic_3_Rt)

figure
%plot(x_plot,ic_1,'g-o',x_plot,ic_2,'b-o',x_plot,ic_3,'c-o',x_plot,ic_4,'m-o',x_plot,ic_5,'r-o')
stem(transpose(x_plot), transpose([ic_1; ic_2; ic_3]),'LineStyle','--')
legend({'$IC_1$','$IC_2$','$IC_3$'},'Interpreter','latex','fontsize',18)
xlabel('Number of data points','FontSize',13)
ylabel('RMSE','FontSize',13)
set(gcf,'Position',[900 400 900 400])
print(gcf,'RMSE_different_data_points_poly_4','-dpdf','-r300','-bestfit')
system(['pdfcrop',' ','RMSE_different_data_points_poly_4','.pdf',' ','RMSE_different_data_points_poly_4','.pdf'])

figure
%plot(x_plot,ic_1_R1,'g-o',x_plot,ic_2_R1,'b-o',x_plot,ic_3_R1,'c-o',x_plot,ic_4_R1,'m-o',x_plot,ic_5_R1,'r-o')
stem(transpose(x_plot), transpose([ic_1_R1; ic_2_R1; ic_3_R1]),'LineStyle','--')
legend({'$IC_1$','$IC_2$','$IC_3$'},'Interpreter','latex','fontsize',18,'Location','southeast')
ylim([min(min(transpose([ic_1_R1; ic_2_R1; ic_3_R1])))*0.99 1])
xlabel('Number of data points','FontSize',13)
ylabel('Correlation $\dot{x}_1$','Interpreter','latex','FontSize',13)
set(gcf,'Position',[900 400 900 400])
print(gcf,'corr_x_different_data_points_poly_4','-dpdf','-r300','-bestfit')
system(['pdfcrop',' ','corr_x_different_data_points_poly_4','.pdf',' ','corr_x_different_data_points_poly_4','.pdf'])

figure
%plot(x_plot,ic_1_R2,'g-o',x_plot,ic_2_R2,'b-o',x_plot,ic_3_R2,'c-o',x_plot,ic_4_R2,'m-o',x_plot,ic_5_R2,'r-o')
stem(transpose(x_plot), transpose([ic_1_R2; ic_2_R2; ic_3_R2]),'LineStyle','--')
legend({'$IC_1$','$IC_2$','$IC_3$'},'Interpreter','latex','fontsize',18,'Location','southeast')
ylim([min(min(transpose([ic_1_R2; ic_2_R2; ic_3_R2])))*0.99 1])
xlabel('Number of data points','FontSize',13)
ylabel('Correlation $\dot{x}_2$','Interpreter','latex','FontSize',13)
set(gcf,'Position',[900 400 900 400])
print(gcf,'corr_y_different_data_points_poly_4','-dpdf','-r300','-bestfit')
system(['pdfcrop',' ','corr_y_different_data_points_poly_4','.pdf',' ','corr_y_different_data_points_poly_4','.pdf'])

figure
%plot(x_plot,ic_1_R3,'g-o',x_plot,ic_2_R3,'b-o',x_plot,ic_3_R3,'c-o',x_plot,ic_4_R3,'m-o',x_plot,ic_5_R3,'r-o')
stem(transpose(x_plot), transpose([ic_1_R3; ic_2_R3; ic_3_R3]),'LineStyle','--')
legend({'$IC_1$','$IC_2$','$IC_3$'},'Interpreter','latex','fontsize',18,'Location','southeast')
ylim([min(min(transpose([ic_1_R3; ic_2_R3; ic_3_R3])))*0.99 1])
xlabel('Number of data points','FontSize',13)
ylabel('Correlation $\dot{x}_3$','Interpreter','latex','FontSize',13)
set(gcf,'Position',[900 400 900 400])
print(gcf,'corr_z_different_data_points_poly_4','-dpdf','-r300','-bestfit')
system(['pdfcrop',' ','corr_z_different_data_points_poly_4','.pdf',' ','corr_z_different_data_points_poly_4','.pdf'])

figure
%plot(x_plot,ic_1_Rt,'g-o',x_plot,ic_2_Rt,'b-o',x_plot,ic_3_Rt,'c-o',x_plot,ic_4_Rt,'m-o',x_plot,ic_5_Rt,'r-o')
stem(transpose(x_plot), transpose([ic_1_Rt; ic_2_Rt; ic_3_Rt]),'LineStyle','--')
legend({'$IC_1$','$IC_2$','$IC_3$'},'Interpreter','latex','fontsize',18,'Location','southeast')
ylim([min(min(transpose([ic_1_Rt; ic_2_Rt; ic_3_Rt])))*0.99 1])
xlabel('Number of data points','FontSize',13)
ylabel('Phase space correlation','FontSize',13)
set(gcf,'Position',[900 400 900 400])
print(gcf,'corr_ph_different_data_points_poly_4','-dpdf','-r300','-bestfit')
system(['pdfcrop',' ','corr_ph_different_data_points_poly_4','.pdf',' ','corr_ph_different_data_points_poly_4','.pdf'])
