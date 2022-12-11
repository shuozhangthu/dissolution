set(0, 'DefaultAxesFontWeight', 'bold', ...
    'DefaultAxesFontSize', 16, ...
    'DefaultAxesFontAngle', 'normal', ... % Not sure the difference here
    'DefaultAxesFontWeight', 'normal', ... % Not sure the difference here
    'DefaultAxesTitleFontWeight', 'normal', ...
    'DefaultAxesTitleFontSizeMultiplier', 1) ;
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesLineWidth', 1.5)
set(0, 'DefaultLineMarkerSize', 6)

clear
close all

load solution_5
figure
set(gcf,'unit','centimeters','position',[1,0,20,15]);
hold on
scatter(1-omega_5,Rmeas_5,60,[0.9412,0.5765,0.1686],'kd','linewidth',1)
scatter(1-omega_5,Rb_solution,60,[0.9216,0.3020,0.2941],'ko','filled','linewidth',1)


load p16
load peterson
load takahashi
omega_peterson=interp1(depth_p16,omega_p16,depth_peterson)
rate_peterson=R_peterson/1000/100*10000/365/24/60/60
scatter(1-omega_peterson,rate_peterson,60,[0.5,0.5,0.5],'o','linewidth',1)


load naviaux_in_situ
scatter(1-omega_naviaux_in_situ,rate_naviaux_in_situ,60,[0.5,0.5,0.5],'+','linewidth',1)

% dim = [.5 .6 .3 .3];
% str = '5 ^oC, D=1\times10^{-22} m^2/s';
% annotation('textbox',dim,'String',str,'fontsize',16,'fontweight','bold','FitBoxToText','on');

load solution_5
scatter(1-omega_5,Rnet_solution,60,[0.4157,0.6902,0.2980],'ks','filled','linewidth',1)



box on
ax = gca;
ax.LineWidth = 1.5;
xlabel('1-\Omega_{calcite}');
ylabel('Rate (mol/m^2/s)');
legend('Naviaux lab 5 ^oC','Model R_{b} 5 ^oC','Peterson {\it in situ}','Naviaux {\it in situ}','Model R_{net} 5 ^oC','location','northwest','fontsize',14);
% set(gca,'XTick',0:1:3);
% set(gca,'YTick',0:1:4);
set(gca,'xscale','log')
set(gca,'yscale','log')

% print('R5_with_peterson.jpeg','-djpeg','-r1200');

% 
% load solution_5
% figure
% set(gcf,'unit','centimeters','position',[1,0,18,15]);
% hold on
% scatter(1-omega_5,Rb_solution,60,[0.9216,0.3020,0.2941],'s','linewidth',2)
% scatter(1-omega_5,Rmeas_5,60,[0.9412,0.5765,0.1686],'d','linewidth',2)
% scatter(1-omega_5,Rnet_solution,60,[0.4157,0.6902,0.2980],'o','linewidth',2)
% xlabel('1-\Omega')
% ylabel('Rate (mol/m^2/s)')
% set(gca,'xscale','log')
% xticks([0.05 0.1 0.5 1])
% xlim([0.05 1])
% set(gca,'yscale','log')
% legend('R_b','R_{meas}','R_{net}','location','southeast')
% box on
% ax = gca;
% ax.LineWidth = 1.5;
% dim = [.2 .6 .3 .3];
% str = '5 ^oC, D=1\times10^{-22} m^2/s';
% annotation('textbox',dim,'String',str,'fontsize',16,'fontweight','bold','FitBoxToText','on');
% % ylim([1e-12 1e-5])
% % print('R5.jpg','-djpeg','-r1200');
% 
% 
% load solution_12
% figure
% set(gcf,'unit','centimeters','position',[1,0,18,15]);
% hold on
% scatter(1-omega_12,Rb_solution,60,[0.9216,0.3020,0.2941],'s','linewidth',2)
% scatter(1-omega_12,Rmeas_12,60,[0.9412,0.5765,0.1686],'d','linewidth',2)
% scatter(1-omega_12,Rnet_solution,60,[0.4157,0.6902,0.2980],'o','linewidth',2)
% xlabel('1-\Omega')
% ylabel('Rate (mol/m^2/s)')
% set(gca,'xscale','log')
% xticks([0.05 0.1 0.5 1])
% xlim([0.05 1])
% set(gca,'yscale','log')
% legend('R_b','R_{meas}','R_{net}','location','southeast')
% box on
% ax = gca;
% ax.LineWidth = 1.5;
% dim = [.2 .6 .3 .3];
% str = '12 ^oC, D=1\times10^{-22} m^2/s, L=7nm';
% annotation('textbox',dim,'String',str,'fontsize',16,'fontweight','bold','FitBoxToText','on');
% % ylim([1e-12 1e-5])
% % print('R12.jpg','-djpeg','-r1200');
% 
% load solution_21
% figure
% set(gcf,'unit','centimeters','position',[1,0,18,15]);
% hold on
% scatter(1-omega_21,Rb_solution,60,[0.9216,0.3020,0.2941],'s','linewidth',2)
% scatter(1-omega_21,Rmeas_21,60,[0.9412,0.5765,0.1686],'d','linewidth',2)
% scatter(1-omega_21,Rnet_solution,60,[0.4157,0.6902,0.2980],'o','linewidth',2)
% xlabel('1-\Omega')
% ylabel('Rate (mol/m^2/s)')
% set(gca,'xscale','log')
% xticks([0.05 0.1 0.5 1])
% xlim([0.05 1])
% set(gca,'yscale','log')
% legend('R_b','R_{meas}','R_{net}','location','southeast')
% box on
% ax = gca;
% ax.LineWidth = 1.5;
% dim = [.2 .6 .3 .3];
% str = '21 ^oC, D=1\times10^{-22} m^2/s, L=7nm';
% annotation('textbox',dim,'String',str,'fontsize',16,'fontweight','bold','FitBoxToText','on');
% % ylim([1e-12 1e-5])
% % print('R21.jpg','-djpeg','-r1200');
% 
% 
% load solution_37
% figure
% set(gcf,'unit','centimeters','position',[1,0,18,15]);
% hold on
% scatter(1-omega_37,Rb_solution,60,[0.9216,0.3020,0.2941],'s','linewidth',2)
% scatter(1-omega_37,Rmeas_37,60,[0.9412,0.5765,0.1686],'d','linewidth',2)
% scatter(1-omega_37,Rnet_solution,60,[0.4157,0.6902,0.2980],'o','linewidth',2)
% xlabel('1-\Omega')
% ylabel('Rate (mol/m^2/s)')
% set(gca,'xscale','log')
% xticks([0.05 0.1 0.5 1])
% xlim([0.05 1])
% set(gca,'yscale','log')
% legend('R_b','R_{meas}','R_{net}','location','southeast')
% box on
% ax = gca;
% ax.LineWidth = 1.5;
% dim = [.2 .6 .3 .3];
% str = '37 ^oC, D=1\times10^{-22} m^2/s, L=7nm';
% annotation('textbox',dim,'String',str,'fontsize',16,'fontweight','bold','FitBoxToText','on');
% % ylim([1e-12 1e-5])
% % print('R37.jpg','-djpeg','-r1200');