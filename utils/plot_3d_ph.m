function plot_3d_ph(xA,vA,tA,xB,vB,tB,x,y,t,xlab,ylab,zlab,gtitle,gname)

    plot3(xA,vA,tA,'LineWidth',1.0,'Color','b')
    hold on
    plot3(xB(1:5:end),vB(1:5:end),tB(1:5:end),'k--','linewidth',3.0)
    hold on
    plot3(x,y,t,'go','MarkerSize',4.0,'MarkerEdgeColor','g','MarkerFaceColor','g')
    hold off
    
    set(gcf,'PaperPositionMode','auto')
    
    title(gtitle,'FontSize',15,'FontName','Helvetica');
    grid on;
    
    set(gca,'Box','on');
    view(-37.5,30)

    labX = xlabel(xlab,'FontSize',12,'FontName','Helvetica');
    labY = ylabel(ylab,'FontSize',12,'FontName','Helvetica');
    labZ = zlabel(zlab,'FontSize',12,'FontName','Helvetica');
    
    %legend({'True phase space','Data-driven dynamics','Data'},'Interpreter','latex','fontsize',12)
    
    print(gcf,gname,'-dpdf','-r300','-bestfit');
    system(['pdfcrop test.pdf test.pdf'])