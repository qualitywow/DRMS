clear; clc; close all;

filename = 'differentSNR_r25_60kHz_b3_K4_ue01.svg';

% please replace the value with the simulation result
% placing order: [JSAC,ZTE,HSCC]
yHSCC =[0,0,0,0.1,0.48,0.73,0.97];
yJSAC =[0,0,0,0,0.16,0.49,0.85];
yZTE =[0,0,0,0,0,0.1,0.45];


x = [-20 -18 -16 -14 -12 -10 -8];
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
%set(gca,'XTick',[1:1:5]);
grid on;
%xticklabels({'-20','-18','-16','-14','-12','-10','-8'});
ylim([0, 1]);
xlabel('SNR','FontSize',12); % , 'Times New Roman','FontSize',12
ylabel('First access rate','FontSize',12); % , 'Times New Roman','FontSize',12
legend('Li [11]', 'Zhang [14]', 'DRMS','FontSize',10, 'Location', 'northwest')


saveas(gcf, filename);