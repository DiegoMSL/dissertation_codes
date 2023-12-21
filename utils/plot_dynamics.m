function plot_dynamics(t_true,x_true,t_iden,x_iden,t_data,x_data,gtitle,xlab,ylab,gname,color,x_order,y_order,xmin,xmax,ymin,ymax)
% function to plot the calculated dynamics together with the identified and 
% the data used to feed the SINDy for a determined time period.
    % t_true = vector [:,1] with the time step used to integrate the original
    % dynamic system.
    % x_true = vector [:,1] with the data calculated for the original dynamics.
    % t_iden = vector [:,1] time step to the identified dynamics.
    % x_iden = vector [:,1] with the data to the identified dynamics.
    % t_data = vector [:,1] with the data time used to feed the SINDy.
    % x_data = vector [:,1] with the data points used to feed the SINDy.
    % gtitle = title of the figure.
    % xlab = label of the x-axis.
    % ylab = label of the y-axis.
    % gname = name of the pdf file created.
    % color = color for the line of the true dynamics plot. (don't use
    % green or black).
    % x_order = scale order of the x-axis.
    % y_order = scale order of the y-axis.
    % x_min = minimum value of the x-axis for the plot.
    % x_max = maximun value of the x-axis for the plot.
    % y_min = minimum value of the y-axis for the plot.
    % y_max = maximun value of the y-axis for the plot.

    fig = figure('Name',gname,'NumberTitle','off');
    plot(t_true,x_true,'LineWidth',4.5,'Color',color) % true dynamics
    hold on
    plot(t_iden(1:10:end),x_iden(1:10:end),'k--','linewidth',2.0) % identified
    hold on
    plot(t_data,x_data,'go','MarkerSize',4.0,'MarkerEdgeColor','g','MarkerFaceColor','g') % samples
    hold off

    %set(gcf,'PaperPositionMode','auto')     % auto scale of the plot
    set(gcf,'Position',[450 200 450 200])
    %set(gcf,'Position',[450 200 225 100]); % manual scale of the plot

    %title(['Component #',num2str(j)]);
    title(gtitle,'FontSize',15,'FontName','Helvetica'); % options for the title

    labX = xlabel(xlab,'FontSize',25,'FontName','Helvetica'); % x-label name options
    labY = ylabel(ylab,'FontSize',25,'FontName','Helvetica'); % y-label name options

    set(gca,'color',[1 1 1]); % background color inside of the box
    set(gcf,'color',[1 1 1]); % background color outside of the box
    set(gca,'FontSize',17.5,'FontName','Helvetica'); % options for the scale numbers

    set(gca,'Box','on');                            % box around graph
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);     % color of the box outline
    set(gca,'TickDir','in','TickLength',[.02 .02]); % tick settings plot box
    set(gca,'XMinorTick','on','YMinorTick','on');   % minor tick settings
    set(gca,'XGrid','off','YGrid','off');           % grid settings

    xlim([xmin xmax]);
    ylim([ymin ymax]);

    % Settings for y axis tick labels and order of magnitude
    %set(gca,'XTickMode','manual','YTickMode','manual')          % preserve tick values for all figure sizes
    %set(gca,'XLimMode','manual','YLimMode','manual')            % preserve axis limits for all figure sizes
    yl = get(gca,'ylim');                                       
    set(gca,'yTick',linspace(yl(1),yl(2),5))         % setting number of tick labels to display          
    BD = y_order;                                    % # of SF before the point in highest tick 
                                                     %label (exception: if highest=1 use 0)
    OM = ceil(log10(yl(2)));                         % ceiling order of magnitude
    ryt=get(gca,'ytick')/10^(OM-BD);                 % redefining tick labels
    % Formating new tick labels
    nyt=cell(size(ryt));
    for i=1:length(ryt)
        nyt{i}=sprintf('% 2.1f',ryt(i));
        % example: '% W.Xf' displays fixed-point notation with X
        % digits after the decimal point, minimum of W characters.
        % The space after the percent inserts a space before the
        % displayed value, giving the same size to + and - numbers. 
    end
    set(gca,'yticklabel',nyt);                                  % setting tick labels
    % Placing order of magnitude
    fs = get(gca,'fontsize');
    set(gca,'units','normalized');
    xl = xlim;
    %text(xl(1),yl(2),sprintf('\\times10^{%d}',OM-BD),'fontsize',fs,'VerticalAlignment','bottom');

    % Settings for x axis tick labels and order of magnitude
    xl = get(gca,'xlim');                                       
    set(gca,'xTick',linspace(xl(1),xl(2),5))     % setting number of tick labels to display          
    BD = x_order;                                % # of SF before the point in highest tick 
                                                 % label (exception: if highest=1 use 0)
    OM = ceil(log10(xl(2)));                     % ceiling order of magnitude
    rxt=get(gca,'xtick')/10^(OM-BD);             % redefining tick labels
    % Formating new tick labels
    nxt=cell(size(rxt));
    for i=1:length(rxt)
        nxt{i}=sprintf('% 2.1f',rxt(i));
        % example: '% W.Xf' displays fixed-point notation with X
        % digits after the decimal point, minimum of W characters.
        % The space after the percent inserts a space before the
        % displayed value, giving the same size to + and - numbers. 
    end
    set(gca,'xticklabel',nxt);                                  % setting tick labels
    % Placing order of magnitude
    fs = get(gca,'fontsize');
    set(gca,'units','normalized');
    xl = xlim;
    %text(xl(2),yl(1),sprintf('\\times10^{%d}',OM-BD),'fontsize',fs,'VerticalAlignment','bottom');
    %text('Units','normalized','VerticalAlignment','bottom','FontSize',fs,'String',...
    %	                        sprintf('\\times10^{%d}',OM-BD),'Position',[0.88 -0.35 0]);
    %legend({'Numerical system','Data driven system','Data used'},'Interpreter','latex','fontsize',10)
    
    print(gcf,gname,'-dpdf','-r300','-bestfit');
    system(['pdfcrop',' ', gname,'.pdf',' ',gname,'.pdf'])