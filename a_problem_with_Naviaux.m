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
set(0, 'DefaultLineMarkerSize', 2)
%%--------------------------------------------------------------------------

load subhas_5;


figure
set(gcf,'Position',[1000, 500, 500,420])
hold on
scatter(1-omega_5,Rmeas_5,40,[0.9216,0.3020,0.2941],'s','MarkerFaceColor',[1,0.4745,0.4745],'linewidth',1)
scatter(1-omega_5,Rmeas_5./(1-omega_5),40,[0.4157,0.6902,0.2980],'^','MarkerFaceColor','(0.7294,0.8627,0.3451)','linewidth',1)
xlabel('1-\Omega','fontsize',12)
ylabel('Rates (mol/m^2/s)','fontsize',12)
set(gca,'xscale','log')
set(gca,'yscale','log')
legend('R_{meas} BET','R_b=R_{meas}/(1-\Omega)','location','northwest')
box on
ax = gca;
ax.LineWidth = 1.5;
% title('5 ^oC');
% print('fig2.jpg','-djpeg','-r1200');



