% 1D transport using MATLAB pdepe                   
clear
set(0, 'DefaultAxesFontWeight', 'bold', ...
    'DefaultAxesFontSize', 18, ...
    'DefaultAxesFontAngle', 'normal', ... % Not sure the difference here
    'DefaultAxesFontWeight', 'bold', ... % Not sure the difference here
    'DefaultAxesTitleFontWeight', 'bold', ...
    'DefaultAxesTitleFontSizeMultiplier', 1) ;
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesLineWidth', 2)
set(0, 'DefaultLineMarkerSize', 6)
%%--------------------------------------------------------------------------
T = 60*60*48;                     % maximum time [s]
L = 1e-7;                     % length [m]

M = 1000;                    % number of timesteps
N = 100;                    % number of nodes  a

t = linspace (0,T,M);    % time discretization
x = linspace (0,L,N);      % space discretization
%

m = 0;
xOde = 0;

Rb=5e-10;
Omega=0.95;
Rnet=Rb*(1-Omega);

Mcc = 26300;                    % mol/m3
D =1e-23;                     % diffusivity [m*m/s]

L=0.3e-9;


pde = @(x,t,u,DuDx) pdeFunc(x,t,u,DuDx,D,Mcc,Rnet);
ic = @(x) icFunc(x);
bc = @(xl,ul,xr,ur,t,v,vdot) bcFunc(xl,ul,xr,ur,t,v,vdot,Rb,Rnet,D, Mcc,L);

ode = @(t,v,vdot,x,u,DuDx) odeFunc(t,v,vdot,x,u,DuDx,D,Mcc,Rb,Rnet,L);

odeic = @() odeIcFunc();

opts.vectorized='off'; % speed up computation
[u,v] = pde1dm(m, pde,ic,bc,x,t,ode, odeic,xOde,opts);


figure
set(gcf,'unit','centimeters','position',[1,0,18,15]);
hold on
load subhas_SIMS.mat;

plot(depth_unreacted, f_unreacted,'k--')
plot(depth_reacted, f_reacted,'k-')

% plot (x*1e9,f13(10,:))        % 
% plot (x*1e9,f13(100,:))        % 
plot ([0 L*1e9 (x(1:end)+L)*1e9],[v(end) v(end) u(end,1:end)],'color',[0.9216,0.3020,0.2941])        % 
% load steady_Dx10
% plot (x*1e9,f13(1000,:))        % 
% set(gca,'xscale','log')
xlim([0 50])
xlabel ('Depth (nm)'); 
ylabel ('f_{13}');
% set(gca,'xscale','log')
% title('D=10^{-22} m^2/s, Rnet=10^{-11} mol/m^2/s, Rb=10^{-9} mol/m^2/s')
box on
ax = gca;
ax.LineWidth = 1.5;


% legend('0.48 hour','4.8 hours','48 hours','steady state','location','best')
legend('Initial, SIMS','48 hours, SIMS','48 hours, modeled','location','southeast')

% print('f13subhas_SIMS.jpg','-djpeg','-r1200');

figure
set(gcf,'unit','centimeters','position',[20,0,18,15]);
plot(t/60/60,v(:))
xlabel('Time (hours)')
ylabel('f_{13}')

box on
ax = gca;
ax.LineWidth = 1.5;

% text(5,0.7,'Subhas et al. (2015), GCA, Fig.5b,\newline21 ^oC, B23-B2\newline\Omega=0.66, R_{meas}=1.67\times10^{-8} mol/m^2/s','FontSize',16)

% figure
% plot(t/60/60,f13(:,1))
% xlabel('Time (hours)')
% ylabel('f^{13}_s (x=0)')
% % print('f13subhas_SIMS_time.jpg','-djpeg','-r1200');



%% functions
function [c,f,s] = pdeFunc(x,t,u,DuDx,d,mcc,rnet)
c = 1;
f = d*DuDx;
s = rnet/mcc*DuDx;
end

function u0 = icFunc(x)
load subhas_SIMS.mat
depth_e=[0;depth_unreacted;1000];
f_e=[f_unreacted(1);f_unreacted;f_unreacted(end)];
u0 = interp1(depth_e/1e9,f_e,x);
end

function [pl,ql,pr,qr] = bcFunc(xl,ul,xr,ur,t,v,vdot,rb,rnet,d, mcc,l)
% pl = -rb/mcc*v(1)+rnet/mcc*ul-l*vdot(1);
% pl = -rb/mcc*(v)+rnet/mcc*(ul);
pl = 2*d*(v-ul)/l;
ql = 1;
pr = ur-1;
qr = 0;
end

function f=odeFunc(t,v,vdot,x,u,DuDx,d,mcc,rb,rnet,l)
f=-rb*v+d*mcc*DuDx+rnet*u-l*mcc*vdot;
end

function v0=odeIcFunc()
v0=0.4516;
end

