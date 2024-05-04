% 1D transport using MATLAB pdepe                   
clear
% close all

set(0, 'DefaultAxesFontWeight', 'normal', ...
    'DefaultAxesFontSize', 14, ...
    'DefaultAxesFontAngle', 'normal', ... % Not sure the difference here
    'DefaultAxesFontWeight', 'normal', ... % Not sure the difference here
    'DefaultAxesTitleFontWeight', 'normal', ...
    'DefaultAxesTitleFontSizeMultiplier', 1) ;
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesLineWidth', 2)
set(0, 'DefaultLineMarkerSize', 6)
%%--------------------------------------------------------------------------
T = 172800;                     % maximum time [s]
L = 1e-6;                     % length [m]

M = 1000;                    % number of timesteps
N = 1000;                    % number of nodes  a

t = linspace (0,T,M);    % time discretization
x = linspace (0,L,N);      % space discretization

load naviaux_5;

for i=1:length(omega_5)
    
Omega=omega_5(i);
Rmeas=Rmeas_5(i);

% Rnet=0;
Mcc = 26300;                    % mol/m3
D = 50e-23;                     % diffusivity [m*m/s]

Rb=[]; Rnet=[];

Rb=[Rmeas:Rmeas:Rmeas/(1-Omega)];
Rnet=Rb*(1-Omega);

f_final=[]; Rm=[];

for j=1:length(Rb)

f13 = pdepe(0,@transfun,@ictransfun,@bctransfun,x,t,odeset,D,Mcc,Rb(j),Rnet(j));
        f_final(j)=(f13(1000,1)+f13(500,1))/2;
Rm(j) = Rb(j)*f_final(j);

end
% 
% 
% figure
% plot(Rb,f_final)
% xlabel('R_b (mol/m^2/s)')
% ylabel('f^{13}_s (x=0)')
% 
% figure
% plot(Rb,Rm)
% hold on
% plot(Rb,Rnet)
% % plot(Rb,Rb)
% plot(Rb,Rb*0+Rmeas)
% xlabel('R_b (mol/m^2/s)')
% ylabel('R_{meas}')
% % set(gca,'xscale','log')
% xlim([1e-10 5e-9])

diff=abs(Rm-Rmeas);

[M(i),I(i)]=min(diff);

Rb_solution(i)=Rb(I(i));
Rnet_solution(i)=Rb_solution(i)*(1-Omega);
f_solution(i)=f_final(I(i));

i

end 

figure
hold on
scatter(1-omega_5(1:end),Rb_solution(1:end),80,'ks','linewidth',1)
scatter(1-omega_5(1:end),Rmeas_5(1:end),80,'ko','linewidth',1)
scatter(1-omega_5(1:end),Rnet_solution(1:end),80,'k','filled','s','linewidth',1)

xlabel('1-\Omega')
ylabel('Rate (mol/m^2/s)')
set(gca,'xscale','log')
set(gca,'yscale','log')
[h,icons] = legend('R_b (model)','R_{meas}','R_{net} (model)','FontSize',16,'Position',[.18 .5 .3 .2])
legend('boxoff')
icons = findobj(icons,'Type','line');
% Find lines that use a marker
icons = findobj(icons,'Marker','none','-xor');
% Resize the marker in the legend
set(icons,'MarkerSize',30);

box on 
ax = gca;
ax.LineWidth = 1.5;
dim = [.2 .6 .3 .3];
str = 'Model solid diffusion\newlineeffects on inferred dissolution rate\newline5 ^oC, 48 hours, D_c=50\times10^{-23} m^2/s';
t=annotation('textbox',dim,'String',str,'FitBoxToText','on');
t.EdgeColor='w';
t.FontSize=16;
t2=annotation('textbox',[0.8 0.1 0.1 0.1],'String','d','FitBoxToText','on');
t2.EdgeColor='w';
t2.FontSize=18;
% t.HorizontalAlignment='center'
xlim([0.01 1])
ylim([1e-12 1e-4])
print('R5_50.jpg','-djpeg','-r1200');

figure
scatter(1-omega_5,f_solution,'s','linewidth',2)
xlabel('1-\Omega')
ylabel('f^{13}_s')
set(gca,'xscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;
title('5 ^oC');
% print('f5.jpg','-djpeg','-r1200');


figure
plot(M)
hold on
plot(Rmeas_5)
set(gca,'yscale','log')

save solution_5

%----------------------functions------------------------------
function [c,f,s] = transfun(x,t,u,dudx,d,mcc,rb,rnet)

c = 1;
f = d*dudx;
s = rnet/mcc*dudx;
end
% --------------------------------------------------------------
function u0 = ictransfun(x,d,mcc,rb,rnet)
% load subhas_full
% depth_e=[0;depth_unreacted;1000];
% f_e=[f_unreacted(1);f_unreacted;f_unreacted(end)];
% u0 = interp1(depth_e/1e9,f_e,x);
u0=1;
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bctransfun(xl,ul,xr,ur,t,d,mcc,rb,rnet)

pl = -rb/mcc*(ul-0.01)+rnet/mcc*ul;
ql = 1;
pr = ur-1;
qr = 0;
end


