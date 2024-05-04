% 1D transport using MATLAB pdepe                   
clear
close all

set(0, 'DefaultAxesFontWeight', 'normal', ...
    'DefaultAxesFontSize', 16, ...
    'DefaultAxesFontAngle', 'normal', ... % Not sure the difference here
    'DefaultAxesFontWeight', 'normal', ... % Not sure the difference here
    'DefaultAxesTitleFontWeight', 'normal', ...
    'DefaultAxesTitleFontSizeMultiplier', 1) ;
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesLineWidth', 2)
set(0, 'DefaultLineMarkerSize', 6)
%%--------------------------------------------------------------------------

load naviaux_5;

Mcc=26300;
v1=Rmeas_5/Mcc;
L=3e-10;
D=1e-22;
Pe1=v1*L/D;

v2=Rmeas_5.*(1-omega_5)/Mcc;
Pe2=v2*L/D;

figure
set(gcf,'unit','centimeters','position',[40,20,18,15]);
hold on
scatter(1-omega_5,Pe1,60,'ks','linewidth',1)
scatter(1-omega_5,Pe2,60,'ko','linewidth',1)
xlabel('1-\Omega','fontsize',18)
ylabel('Peclet number','fontsize',18)
set(gca,'xscale','log')
set(gca,'yscale','log')
legend('v_{net}= R_{meas}/M_{cc}','v_{net}= R_{meas}(1-\Omega)/M_{cc}','location','northwest')
box on
ax = gca;
ax.LineWidth = 1.5;
% title('5 ^oC');
print('peclet.jpg','-djpeg','-r1200');
% 


