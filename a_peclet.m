% 1D transport using MATLAB pdepe                   
clear
close all

set(0, 'DefaultAxesFontWeight', 'bold', ...
    'DefaultAxesFontSize', 12, ...
    'DefaultAxesFontAngle', 'normal', ... % Not sure the difference here
    'DefaultAxesFontWeight', 'normal', ... % Not sure the difference here
    'DefaultAxesTitleFontWeight', 'normal', ...
    'DefaultAxesTitleFontSizeMultiplier', 1) ;
set(0, 'DefaultLineLineWidth', 1);
set(0, 'DefaultAxesLineWidth', 1)
set(0, 'DefaultLineMarkerSize', 6)
%%--------------------------------------------------------------------------

load subhas_5;

Mcc=26300;
v1=Rmeas_5/Mcc;
L=3e-9;
D=5e-23;
Pe1=v1*L/D;

v2=Rmeas_5.*(1-omega_5)/Mcc;
Pe2=v2*L/D;

figure
set(gcf,'Position',[1000, 500, 500,420])
hold on
scatter(1-omega_5,Pe1,40,[0.9216,0.3020,0.2941],'s','MarkerFaceColor',[1,0.4745,0.4745],'linewidth',1)
scatter(1-omega_5,Pe2,40,[0.4157,0.6902,0.2980],'^','MarkerFaceColor','(0.7294,0.8627,0.3451)','linewidth',1)
xlabel('1-\Omega','fontsize',12)
ylabel('Peclet number','fontsize',12)
set(gca,'xscale','log')
set(gca,'yscale','log')
legend('Pe (R_{meas})','Pe (R_{meas}(1-\Omega))','location','northwest')
box on
ax = gca;
ax.LineWidth = 1.5;
% title('5 ^oC');
% print('peclet.jpg','-djpeg','-r1200');
% 


