clear; clc; close all;

% output filename
filename = 'differentEK_snr-8_r25_60kHz_K8_ue01.svg';

% please replace the value with the simulation result
% placing order: [JSAC,ZTE,HSCC]
bSingle = [0.84,0.87,0.96]; % |{e/K}| = 1
bMultiple = [0.77,0.79,0.93]; % |{e/K}| = 2
b = [bSingle; bMultiple]


figure();
% color of lines
co = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.4660, 0.6740, 0.1880]];
p = bar(b,0.7);

for i = 1:size(b,2)
    p(i).FaceColor = co(i,:);
end

xticklabels({'[\epsilon/K] = \{0\}', '[\epsilon/K] = \{0, 1\}'});
ylim([0, 1]);
xlabel('Number of users','FontSize',1); % , 'Times New Roman','FontSize',12
ylabel('First access rate','FontSize',12); % , 'Times New Roman','FontSize',12
lgd = legend('Li [11]', 'Zhang [14]', 'DRMS','FontSize',10,'Location', 'northoutside')
lgd.NumColumns = 3;

saveas(gcf, filename);

