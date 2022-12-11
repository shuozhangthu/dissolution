% 1D transport using MATLAB pdepe
clear
% close all

set(0, 'DefaultAxesFontWeight', 'bold', ...
    'DefaultAxesFontSize', 14, ...
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

load subhas_5;

m = 0;
xOde = 0;
Mcc = 26300;                    % mol/m3
D = 10e-23;                     % diffusivity [m*m/s]

L=3e-9;

for i=1:13

    Omega=omega_5(i);
    Rmeas=Rmeas_5(i);


    Rb=[]; Rnet=[];

    if i==7 || i==9
        Rb=[Rmeas:0.1*Rmeas:Rmeas/(1-Omega)];
    elseif i==8
        Rb=[Rmeas:0.1*Rmeas:Rmeas/(1-Omega)];
    elseif i==13 
        Rb=[Rmeas:0.1*Rmeas:Rmeas/(1-Omega)];
    else
        Rb=[Rmeas:Rmeas:Rmeas/(1-Omega)];
    end


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

        f_final(j)=v(end);

        Rm(j) = Rb(j)*f_final(j);

    end


    diff=abs(Rm-Rmeas);

    [M(i),I(i)]=min(diff);

    Rb_solution(i)=Rb(I(i));
    Rnet_solution(i)=Rb_solution(i)*(1-Omega);
    f_solution(i)=f_final(I(i));

    i

end

for i=14:length(omega_5)

    Omega=omega_5(i);
    Rmeas=Rmeas_5(i);

    Rb=[]; Rnet=[];


    Rb=[Rmeas:Rmeas:Rmeas/(1-Omega)];


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

        f_final(j)=v(end);

        Rm(j) = Rb(j)*f_final(j);

    end


    diff=abs(Rm-Rmeas);

    [M(i),I(i)]=min(diff);

    Rb_solution(i)=Rb(I(i));
    Rnet_solution(i)=Rb_solution(i)*(1-Omega);
    f_solution(i)=f_final(I(i));

    i

end

figure
hold on
scatter(1-omega_5,Rb_solution,'s','linewidth',1.5)
scatter(1-omega_5,Rmeas_5,'o','linewidth',1.5)
scatter(1-omega_5,Rnet_solution,'d','linewidth',1.5)
xlabel('1-\Omega')
ylabel('Rates (mol/m^2/s)')
set(gca,'xscale','log')
set(gca,'yscale','log')
legend('R_b','R_{meas}','R_{net}','location','southeast')
box on
ax = gca;
ax.LineWidth = 1.5;
dim = [.2 .6 .3 .3];
str = '5 ^oC, D=5\times10^{-23} m^2/s, L=3 nm';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
xlim([0.01 1])
ylim([1e-12 1e-4])
% print('R5_D5.jpg','-djpeg','-r1200');

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



