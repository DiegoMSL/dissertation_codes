% Rosller example

clear; close all; clc

%% Rosller attractor
% x' = -y-z
% y' = x + alpha*y
% z' = beta + z*(x-gamma)
rng(16081995)
alpha = 0.1;
beta = 0.1;
gamma = 14;
param = [alpha; beta; gamma];
x0 = [-8; 8; 0];

equation = @(t,x)rossler(t,x,param);

%% Parameters SINDy

lambda = 0.05;
polyorder = 5;
usesine = 0;
n = 3;
dt_jump = 100;         % jump in the int_dt for the sparse dataset


% Integrate
tspan=[.001:.001:25];
N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,x]=ode45(equation,tspan,x0,options);

% Sparse Dataset
cont = 1;
for i=1:dt_jump:(N)
    data_t(cont,1) = t(i,1);
    data_x(cont,:) = x(i,:);
    cont = cont+1;
end

t = data_t;
x = data_x;

%% compute Derivative
eps = 0.5;
for i=1:length(x)
    dx(i,:) = rossler(0,x(i,:),param);
end
dx = dx + eps*randn(size(dx));

%% pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(x,n,polyorder,usesine);
m = size(Theta,2);

%% compute Sparse regression: sequential least squares

Xi = sparsifyDynamics(Theta,dx,lambda,n);
poolDataLIST({'x','y','z'},Xi,n,polyorder,usesine);

%% Rosller for T\in[0,50]
tspan = [0 100];
[tA,xA]=ode45(equation,tspan,x0,options);   % true model
[tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);  % approximate

figure
subplot(1,2,1)
title('True attractor')
dtA = [0; diff(tA)];
color_line3(xA(:,1),xA(:,2),xA(:,3),dtA,'LineWidth',1.5);
hold on;
plot3(data_x(:,1),data_x(:,2),data_x(:,3),'ro','MarkerSize',3.0,'MarkerEdgeColor','r','MarkerFaceColor','r')%,'DisplayName','Data')
%legend('','Data')
view(27,16)
set(gca,'FontSize',13.0,'FontName','Helvetica');
grid on
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)
zlabel('z','FontSize',13)
xlim([-25 25]);
ylim([-25 25]);
zlim([0 40]);
subplot(1,2,2)
title('Reconstructed attractor')
dtB = [0; diff(tB)];
color_line3(xB(:,1),xB(:,2),xB(:,3),dtB,'LineWidth',1.5);
view(27,16)
set(gca,'FontSize',13.0,'FontName','Helvetica');
grid on
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)
zlabel('z','FontSize',13)
xlim([-25 25]);
ylim([-25 25]);
zlim([0 40]);
print(gcf,'fig1','-dpdf','-r600','-bestfit');
system(['pdfcrop',' ','fig1','.pdf',' ','fig1','.pdf'])

% Rosller for t=50, dynamo view
%figure
%plot(tA,xA(:,1),'r','LineWidth',2.0), hold on
%plot(tB,xB(:,1),'k--','LineWidth',1.5)
%plot(data_t,data_x(:,1),'go','MarkerSize',3.0,'MarkerEdgeColor','g','MarkerFaceColor','g')
%grid on
%ylabel('x','FontSize',13)
%set(gca,'FontSize',13)
%print(gcf,'fig2','-dpdf','-r600','-bestfit');
%system(['pdfcrop',' ','fig2','.pdf',' ','fig2','.pdf'])
plot_dynamics(tA,xA(:,1),tB,xB(:,1),data_t,data_x(:,1),'','time','x','time_series_x','r',2,1,0,100,-25,25)

plot_dynamics(tA,xA(:,2),tB,xB(:,2),data_t,data_x(:,2),'','time','y','time_series_y','b',2,1,0,100,-25,25)

plot_dynamics(tA,xA(:,3),tB,xB(:,3),data_t,data_x(:,3),'','time','z','time_series_z','c',2,1,0,100,0,40)

%plot(tA,xA(:,2),'r','LineWidth',1.5), hold on
%plot(tB,xB(:,2),'k--','LineWidth',1.5)
%plot(data_t,data_x(:,2),'go','MarkerSize',3.0,'MarkerEdgeColor','g','MarkerFaceColor','g')
%grid on
%xlabel('time','FontSize',13)
%ylabel('y','FontSize',13)


%plot(tA,xA(:,3),'r','LineWidth',1.5), hold on
%plot(tB,xB(:,3),'k--','LineWidth',1.5)
%plot(data_t,data_x(:,3),'go','MarkerSize',3.0,'MarkerEdgeColor','g','MarkerFaceColor','g')
%grid on
%ylabel('z','FontSize',13)

%set(gcf,'PaperPositionMode','auto')     % auto scale of the plot
%set(gcf,'Position',[900 500 900 500]); % manual scale of the plot
%print(gcf,'fig2','-dpdf','-r600','-bestfit');
%system(['pdfcrop',' ','fig2','.pdf',' ','fig2','.pdf'])

%% Rosller for T\in[0,500]
%tspan = [0 500];
%options = odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,n));
%[tA,xA]=ode45(equation,tspan,x0,options);   % true model
%[tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0,options);  % approximate

%figure
%dtA = [0; diff(tA)];
%color_line3(xA(:,1),xA(:,2),xA(:,3),dtA,'LineWidth',1.5);
%view(27,16)
%grid off
%xlabel('x','FontSize',13)
%ylabel('y','FontSize',13)
%zlabel('z','FontSize',13)
%set(gcf,'PaperPositionMode','auto');
%set(gca,'Box','off');                            % box around graph
%set(gca,'TickDir','in','TickLength',[.000000002 .0000000000002]);      % tick settings
%set(gca,'FontSize',0.00001,'FontName','Helvetica');
%set(gca,'XColor', 'none','YColor','none','ZColor','none')
%gname=['test'];
%print(gcf,gname,'-dpdf','-r600','-bestfit');
%system(['pdfcrop',' ',gname,'.pdf',' ',gname,'.pdf'])