close all
font_size = 15;
line_size = 15;
line_width = 2;

Re = cst.Re;

figure
plot3(SC_r_g_initial_hat(:,1),SC_r_g_initial_hat(:,2),SC_r_g_initial_hat(:,3),'rx','Linewidth',line_width);
hold on
plot3(SC_r_g_hat(:,1),SC_r_g_hat(:,2),SC_r_g_hat(:,3),'mx','Linewidth',line_width);
hold on
grid on
[F_g1, F_g2, F_g3] = sphere(36);
colormap(white)
surf(F_g1*Re, F_g2*Re, F_g3*Re)
load coast.mat
[x_coast, y_coast, z_coast] = sph2cart(long*(pi/180), lat*(pi/180), Re);
plot3(x_coast, y_coast, z_coast, 'k')
axis('equal')
% title('Spacecraft Orbit')
xlabel('$g^1$ (m)','fontsize',font_size,'Interpreter','latex');
ylabel('$g^2$ (m)','fontsize',font_size,'Interpreter','latex');
zlabel('$g^3$ (m)','fontsize',font_size,'Interpreter','latex');
view([1 1 1])
print('plots/orbit','-depsc','-r720');

figure
plot(t/T,bias_hat,'b-','Linewidth',line_width);
hold on
% title('Absolute Error vs. $h$','fontsize',font_size,'Interpreter','latex');
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\hat{b}$ (m)','fontsize',font_size,'Interpreter','latex');
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on
exportfig(gcf,'plots/bias_hat.eps','width',11,'Height',8.5,'fontmode','fixed','fontsize',18,'Color','cmyk','LineWidth',line_width);
