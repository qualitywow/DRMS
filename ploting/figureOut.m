clear; clc; close all;

filename = 'differentUE_b6_snr-8_r25_60kHz_K8.svg';

% please replace the value with the simulation result
% placing order: [JSAC,ZTE,HSCC]
yHSCC = [0.885,0.8561,0.6667,0.1944,0.0375];
yJSAC = [0.76,0.6925,0.5694,0.2191,0.0448];
yZTE = [0.795,0.6425,0.5102,0.1637,0.0262];



x = [1 2 3 4 5];
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
set(gca,'XTick',[1:1:5]);
grid on;
xticklabels({'2','4','8','16','32'});
ylim([0, 1]);
xlabel('Number of users','FontSize',12); % , 'Times New Roman','FontSize',12
ylabel('First access rate','FontSize',12); % , 'Times New Roman','FontSize',12
legend('Li [11]', 'Zhang [14]', 'DRMS','FontSize',10,'Location', 'northeast')


saveas(gcf, filename);
