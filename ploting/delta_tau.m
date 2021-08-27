clear;clc;close all;


load 'PPP.mat'

a0 = [188];
a = [388 449 93];

eK = [0];
% eK = [0, 29];
veK = ones(size(eK)) * 533.4;


pdp = PPP{1};

% plot(0:length(pdp)-1,pdp,'s-');
figure();

% u1 other peaks
stem(a,pdp(a),'LineStyle','--',...
     'LineWidth', 1.2, ...
     'MarkerSize', 8, 'Color', [0 0.4470 0.7410]);
hold on;

% u1 first peak
stem(a0,pdp(a0),'LineStyle','--',...
     'LineWidth', 1.2, ...
     'MarkerSize', 8, ...
     'Color', [0 0.4470 0.7410], ...
     'MarkerFaceColor', [0 0.4470 0.7410]);
hold on;

% [0.9290, 0.6940, 0.1250] yellow

% plot delta_1, delta_1', adjust by array eK
stem(eK, veK,'LineStyle','-',...
     'LineWidth', 1.5, 'Marker', '^', ...
     'MarkerSize', 8, 'Color', [0 0.4470 0.7410])
hold on;

% plot u2 peak
% stem(377,(pdp(a0)-50),'LineStyle','--',...
%      'LineWidth', 1.2, ...
%      'MarkerSize', 8, ...
%      'Color', [0.9290, 0.6940, 0.1250], ...
%      'MarkerFaceColor', [0.9290, 0.6940, 0.1250]);
% hold on;

% plot delta_2
% stem([0,220], ones(size([0,220]))*(pdp(a0)-50),'LineStyle','-',...
%      'LineWidth', 1.5, 'Marker', '^', ...
%      'MarkerSize', 8, 'Color', [0.9290, 0.6940, 0.1250])
%  hold off;

% legend for u1_tau1, u1_tau2
l = legend('other peaks of the peak group',...
        'first peak of the peak group',...
        'candidate delay index', 'FontSize', 11, 'Location', 'none');

% legend for u1_u2
% l = legend('first peak of the peak group u1',...
%            'candidate delay indices of u1',...
%            'first peak of the peak group u2', ...
%            'candidate delay indices of u2',...
%            'FontSize', 11, 'Location', 'none', 'Interpreter', 'tex');
lgd.NumColumns = 2;

l.Position = [0.57  0.27 0.1778 0.1957];
xlim([-50 600])
ylim([0 600])
xticks(0:100:600)
xlabel('Delay index','FontSize',12,'Interpreter','none');%'Position', [420 -0.08],
ylabel('Correlation','FontSize',12,'Interpreter','none');%'Position', [-200 0.5],