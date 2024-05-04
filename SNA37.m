% 1D transport using MATLAB pdepe
clear
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
T = 172800;                     % maximum time [s]
L = 1e-6;                     % length [m]

M = 1000;                    % number of timesteps
N = 1000;                    % number of nodes  a

t = linspace (0,T,M);    % time discretization
x = linspace (0,L,N);      % space discretization

load naviaux_37;

for i=1:length(omega_37)
    
    Omega=omega_37(i);
    Rmeas=Rmeas_37(i);
    
    % Rnet=0;
    Mcc = 26300;                    % mol/m3
    D =10e-23;                     % diffusivity [m*m/s]
    
    Rb=[]; Rnet=[];
    
    Rb=[Rmeas:Rmeas:Rmeas/(1-Omega)];
    Rnet=Rb*(1-Omega);
    
    f_final=[]; Rm=[];

    for j=1:length(Rb)
           
        f13 = pdepe(0,@transfun,@ictransfun,@bctransfun,x,t,odeset,D,Mcc,Rb(j),Rnet(j));
        f_final(j)=(f13(1000,1)+f13(500,1))/2;
%         f_final(j)=f13(1000,1);
        Rm(j) = Rb(j)*f_final(j);
       
    end
    
  
    diff=abs(Rm-Rmeas);
    
    [M(i),I]=min(diff);
    
    Rb_solution(i)=Rb(I);
    Rnet_solution(i)=Rb_solution(i)*(1-Omega);
    f_solution(i)=f_final(I);
    
    i
end

figure
hold on
scatter(1-omega_37,Rb_solution,'s','linewidth',2)
scatter(1-omega_37,Rmeas_37,'o','linewidth',2)
scatter(1-omega_37,Rnet_solution,'d','linewidth',2)
xlabel('1-\Omega')
ylabel('Rates (mol/m^2/s)')
set(gca,'xscale','log')
set(gca,'yscale','log')
legend('R_b','R_{meas}','R_{net}','location','southeast')
box on
ax = gca;
ax.LineWidth = 1.5;
dim = [.2 .6 .3 .3];
str = '37 ^oC, D=10\times10^{-23} m^2/s';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
xlim([0.01 1])
ylim([1e-12 1e-4])
print('R37_10.jpg','-djpeg','-r1200');

figure
scatter(1-omega_37,f_solution,60,[0.9216,0.3020,0.2941],'s','linewidth',2)
xlabel('1-\Omega')
ylabel('f^{13}_s')
set(gca,'xscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;
str = '37 ^oC, D=2.5\times10^{-23} m^2/s';
% print('f37.jpg','-djpeg','-r1200');

figure
plot(M)
hold on
plot(Rmeas_37)
set(gca,'yscale','log')

save solution_37

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


