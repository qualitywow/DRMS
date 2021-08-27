clear; clc; close all;

filename = 'differentE_ue_01_snr-8_r25_60kHz_K4.svg';

x = [0,0.5,1.0,1.5,2];

yHSCC = [0.99,0.99,0.98,0.98,0.98];
yJSAC = [0.86,0.91,0.9,0.87,0.9];
yZTE = [0.54,0.45,0.62,0.47,0.39];


figure();
co = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.4660, 0.6740, 0.1880]];
plot(x, yJSAC,'s-',...
    'Color', co(1,:),'MarkerFaceColor', co(1,:), ...
    'MarkerSize', 8); hold on;
plot(x, yZTE,'o-',...
    'Color', co(2,:),'MarkerFaceColor', co(2,:), ...
    'MarkerSize', 7); hold on;
plot(x, yHSCC,'^-',...
    'Color', co(3,:),'MarkerFaceColor', co(3,:), ...
    'MarkerSize', 7); hold on;
set(gca,'XTick',[0:0.5:2]);
grid on;

xlim([0,2]);
%xticklabels({'0','0.5','1.0','1.5','2'});

ylim([0, 1]);
xlabel('Normalized CFO \epsilon','FontSize',12); % , 'Times New Roman','FontSize',12
ylabel('First access rate','FontSize',12); % , 'Times New Roman','FontSize',12
legend('Li [11]', 'Zhang [14]', 'DRMS','FontSize',10, 'Location', 'southeast');


saveas(gcf, filename);