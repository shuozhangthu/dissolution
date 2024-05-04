% 1D transport using MATLAB pdepe                   
clear
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
T = 60*60*48;                     % maximum time [s]
L = 1e-6;                     % length [m]

M = 1000;                    % number of timesteps
N = 1000;                    % number of nodes  a

t = linspace (0,T,M);    % time discretization
x = linspace (0,L,N);      % space discretization

Rb=4.94e-10;
Omega=0.95;
Rnet=Rb*(1-Omega);
% Rnet=0;
Mcc = 26300;                    % mol/m3
D = 1e-23;                     % diffusivity [m*m/s]


f13 = pdepe(0,@transfun,@ictransfun,@bctransfun,x,t,odeset,D,Mcc,Rb,Rnet);

figure
set(gcf,'unit','centimeters','position',[1,0,18,15]);
hold on
load subhas_full
plot(depth_unreacted, f_unreacted,'k--')
plot(depth_reacted, f_reacted,'k-')

% plot (x*1e9,f13(10,:))        % 
% plot (x*1e9,f13(100,:))        % 
plot (x*1e9,f13(1000,:),'color',[0.9216,0.3020,0.2941])        % 
% load steady_Dx10
% plot (x*1e9,f13(1000,:))        % 
% set(gca,'xscale','log')
xlim([0 50])
xlabel ('Depth (nm)'); 
ylabel ('f_{ic}');
% set(gca,'xscale','log')
% title('D=10^{-22} m^2/s, Rnet=10^{-11} mol/m^2/s, Rb=10^{-9} mol/m^2/s')
box on
ax = gca;
ax.LineWidth = 1.5;


% legend('0.48 hour','4.8 hours','48 hours','steady state','location','best')
legend('Initial, SIMS','48 hours, SIMS','48 hours, modeled','FontSize',16,'Position',[.37 .2 .3 .2])
legend('boxoff')


dim = [.35 .3 .3 .3];
str = 'Isotope ratio in solid calcite\newline21 ^oC, 48 hours, D_c=1\times10^{-23} m^2/s\newline\Omega=0.95, R_b= 5\times 10^{-10} mol/m^2/s';
t1=annotation('textbox',dim,'String',str,'FitBoxToText','on');
t1.EdgeColor='w';
t1.FontSize=16;
% t2=annotation('textbox',[0.8 0.1 0.1 0.1],'String','d','FitBoxToText','on');
% t2.EdgeColor='w';
% t2.FontSize=18;

print('f13subhas_SIMS.jpg','-djpeg','-r1200');

figure
set(gcf,'unit','centimeters','position',[1,0,18,15]);
plot(t/60/60,f13(:,1))
xlabel('Time (hours)')
ylabel('f_{13}')

box on
ax = gca;
ax.LineWidth = 1.5;

text(5,0.7,'Subhas et al. (2015), GCA, Fig.5b,\newline21 ^oC, B23-B2\newline\Omega=0.66, R_{meas}=1.67\times10^{-8} mol/m^2/s','FontSize',16)

% figure
% plot(t/60/60,f13(:,1))
% xlabel('Time (hours)')
% ylabel('f^{13}_s (x=0)')
% % print('f13subhas_SIMS_time.jpg','-djpeg','-r1200');


%----------------------functions------------------------------
function [c,f,s] = transfun(x,t,u,dudx,d,mcc,rb,rnet)

c = 1;
f = d*dudx;
s = rnet/mcc*dudx;
end
% --------------------------------------------------------------
function u0 = ictransfun(x,d,mcc,rb,rnet)
load subhas_full
depth_e=[0;depth_unreacted;1000];
f_e=[f_unreacted(1);f_unreacted;f_unreacted(end)];
u0 = interp1(depth_e/1e9,f_e,x);
% u0=1;
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bctransfun(xl,ul,xr,ur,t,d,mcc,rb,rnet)

pl = -rb/mcc*(ul-0.01)+rnet/mcc*ul;
ql = 1;
pr = ur-1;
qr = 0;
end


