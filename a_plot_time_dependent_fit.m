clear
set(0, 'DefaultAxesFontWeight', 'normal', ...
    'DefaultAxesFontSize', 12, ...
    'DefaultAxesFontAngle', 'normal', ... % Not sure the difference here
    'DefaultAxesFontWeight', 'normal', ... % Not sure the difference here
    'DefaultAxesTitleFontWeight', 'normal', ...
    'DefaultAxesTitleFontSizeMultiplier', 1) ;
set(0, 'DefaultLineLineWidth', 1);
set(0, 'DefaultAxesLineWidth', 1)
set(0, 'DefaultLineMarkerSize', 6)


figure
hold on
set(gcf,'Position',[1000, 500, 500,420])
xlabel('Time (hours)')
ylabel('% solid dissolved')
box on
ax = gca;
ax.LineWidth = 1.5;

load percent_dissolved_fast.mat
scatter(time,percent_dissolved,'s','linewidth',2)
load percent_dissolved_slow.mat
scatter(time,percent_dissolved,'o','linewidth',2)

set(gca,'colororderindex',1)

a_fit_time_dependent_data_fast;
plot(t/60/60,cumsum(Rb_solution*v(:)*BET*100*360)*100+0.1,'-')

a_fit_time_dependent_data_slow;
plot(t/60/60,cumsum(Rb_solution*v(:)*BET*100*360)*100,'-')

legend('Experiment, 20-53 \mum','Experiment, 70-100 \mum','Model, 20-53 \mum','Model, 70-100 \mum','location','northwest');

% print('percent_dissolved_fit.jpg','-djpeg','-r1200');
