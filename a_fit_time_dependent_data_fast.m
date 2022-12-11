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

%% experimental data and parameters

load percent_dissolved_fast.mat

time_selected=time(3:end);
percent_dis_selected=percent_dissolved(3:end);

rate=0.0023; %g/g/day

BET=0.1520;%m2/g

Rmeas=rate/BET/100/24/60/60;%mol/m2/s

% Rmeas=Rmeas/10;

T = 60*60*100;                     % maximum time [s]
L = 1e-7;                     % length [m]

M = 1000;                    % number of timesteps
N = 100;                    % number of nodes  a

t = linspace (0,T,M);    % time discretization
x = linspace (0,L,N);      % space discretization
%

m = 0;
xOde = 0;


Omega=0.83;

Mcc = 26300;                    % mol/m3
D =1e-22;                     % diffusivity [m*m/s]
D=0;
L=3e-9;


%% find solution
Rb=[]; Rnet=[];

% Rb=[Rmeas:0.5*Rmeas:Rmeas/(1-Omega)*2];

Rb=Rmeas/(1-Omega);

Rnet=Rb*(1-Omega);

f_final=[]; Rm=[];

for j=1:length(Rb)

    pde = @(x,t,u,DuDx) pdeFunc(x,t,u,DuDx,D,Mcc,Rnet(j));
    ic = @(x) icFunc(x);
    bc = @(xl,ul,xr,ur,t,v,vdot) bcFunc(xl,ul,xr,ur,t,v,vdot,Rb(j),Rnet(j),D,Mcc,L);

    ode = @(t,v,vdot,x,u,DuDx) odeFunc(t,v,vdot,x,u,DuDx,D,Mcc,Rb(j),Rnet(j),L);

    odeic = @() odeIcFunc();

    opts.vectorized='off'; % speed up computation
    [u,v] = pde1dm(m, pde,ic,bc,x,t,ode, odeic,xOde,opts);
    %         f_final(j)=(f13(1000,1)+f13(500,1))/2;
    f_final(j)=v(end);
    Rm(j) = Rb(j)*f_final(j);
    j
end


diff=abs(Rm-Rmeas);

[M,I]=min(diff);

Rb_solution=Rb(I)
Rnet_solution=Rb_solution*(1-Omega);

f_solution=f_final(I);


%% plot solution
pde = @(x,t,u,DuDx) pdeFunc(x,t,u,DuDx,D,Mcc,Rnet_solution);
ic = @(x) icFunc(x);
bc = @(xl,ul,xr,ur,t,v,vdot) bcFunc(xl,ul,xr,ur,t,v,vdot,Rb_solution, Rnet_solution,D,Mcc,L);

ode = @(t,v,vdot,x,u,DuDx) odeFunc(t,v,vdot,x,u,DuDx,D,Mcc,Rb_solution,Rnet_solution,L);

odeic = @() odeIcFunc();

opts.vectorized='off'; % speed up computation

[u,v] = pde1dm(m, pde,ic,bc,x,t,ode, odeic,xOde,opts);

% 
% figure
% set(gcf,'unit','centimeters','position',[1,0,18,15]);
% hold on
% plot(t/60/60,v(:))
% % scatter(t/60/60,u(:,1))
% xlabel('Time (hours)')
% ylabel('f_{s}')
% box on
% ax = gca;
% ax.LineWidth = 1.5;
% text(10,0.8,'Naviaux et al. (2019), GCA, Fig.1a,\newline21 ^oC, 20-53 \mum size\newline\Omega=0.83, R_{b}=9.7\times10^{-9}, R_{net}=1.6\times10^{-9} mol/m^2/s\newline\newlineD=0 m^2/s, L= 3nm','FontSize',16)
% 
% % % % print('fs_time_fast.jpg','-djpeg','-r1200');
% 
% figure;
% plot([0 L*1e9 (x+L)*1e9], [v(end) v(end) u(end,:)],'color',[0.9216,0.3020,0.2941]);
% % plot(x*1e9, u(end,:),'color',[0.9216,0.3020,0.2941]);
% set(gcf,'unit','centimeters','position',[20,0,18,15]);
% xlabel('Depth (nm)');
% xlim([0 50])
% ylabel('f_{s} and f_{b} at 100 hours');
% box on
% ax = gca;
% ax.LineWidth = 1.5;
% text(10,0.8,'Naviaux et al. (2019), GCA, Fig.1a,\newline21 ^oC, 20-53 \mum size\newline\Omega=0.83, R_{b}=9.7\times10^{-9}, R_{net}=1.6\times10^{-9} mol/m^2/s\newline\newlineD=0 m^2/s, L= 3nm','FontSize',16)
% 
% % print('fb_depth_fast.jpg','-djpeg','-r1200');
% 
% figure
% set(gcf,'unit','centimeters','position',[40,0,18,15]);
% plot(t/60/60,cumsum(Rb_solution*v(:)*BET*100*360)*100+0.1)
% hold on
% scatter(time,percent_dissolved,'linewidth',2)
% legend('model','Experimental data','location','southeast');
% xlabel('Time (hours)')
% ylabel('% dissolved')
% box on
% ax = gca;
% ax.LineWidth = 1.5;
% text(10,0.8,'Naviaux et al. (2019), GCA, Fig.1a,\newline21 ^oC, 20-53 \mum size\newline\Omega=0.83, R_{b}=9.7\times10^{-9}, R_{net}=1.6\times10^{-9} mol/m^2/s\newline\newlineD=0 m^2/s, L= 3nm','FontSize',16)

% print('percent_dissolved_fast.jpg','-djpeg','-r1200');
% 
% figure
% plot(diff)
% title('diff')

%% functions
function [c,f,s] = pdeFunc(x,t,u,DuDx,d,mcc,rnet)
c = 1;
f = d*DuDx;
s = rnet/mcc*DuDx;
end

function u0 = icFunc(x)
u0 = 1;
end

function [pl,ql,pr,qr] = bcFunc(xl,ul,xr,ur,t,v,vdot,rb,rnet,d,mcc,l)
% pl = -rb/mcc*v(1)+rnet/mcc*ul-l*vdot(1);
% pl = -rb/mcc*ul+rnet/mcc*ul;
pl = 2*d*(v-ul)/l;
ql = 1;
pr = ur-1;
qr = 0;
end

function f=odeFunc(t,v,vdot,x,u,DuDx,d,mcc,rb,rnet,l)
f=-rb*v(1)+d*mcc*DuDx+rnet*u-l*mcc*vdot(1);
end

function v0=odeIcFunc()
v0=1;
end




