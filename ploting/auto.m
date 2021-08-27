% plot autocorrelation of Zadoff-Chu sequence

clc; clear; close all;
Nzc = 839; u = 1;
Ncs = 173; K = 8;

% color for line and marker
% blue for dark and yellow for light
co = {[0, 0.4470, 0.7410], [0.9290, 0.6940, 0.1250]};

x1 = zadoffChuSeq(1, Nzc);
x2 = circshift(x1, -Ncs);

for m = 0:Nzc-1
    PDP1(m+1) = abs(sum(x1 .* conj(circshift(x1, -m)))/Nzc);
    PDP2(m+1) = abs(sum(x2 .* conj(circshift(x1, -m)))/Nzc);
end

figure();

% 'Color',co{1} change line color
% 'MarkerFaceColor',co{1} change marker color
p1 = plot(0:Nzc-1, PDP1,'o-','LineWidth',1.2,'MarkerSize',8,'Color',co{1}); hold on;
p2 = plot(0:Nzc-1, PDP2,'^-','LineWidth',1.2,'MarkerSize',8,'Color',co{2}); hold on;

% Show marker per 100 data points
p1.MarkerIndices =  1:100:length(PDP1); 
p2.MarkerIndices =  1:100:length(PDP2);

% Format Y for data tips
p1.DataTipTemplate.DataTipRows(end).Format = '%.4f';
p2.DataTipTemplate.DataTipRows(end).Format = '%.4f';

% Specify some data tips
dt_u1_1 = datatip(p1,'DataIndex',1,'FontSize',12);
dt_u1_2 = datatip(p1,'DataIndex',301,'FontSize',12);
dt_u2_1 = datatip(p2,'DataIndex',Ncs+1,'FontSize',12);
dt_u2_2 = datatip(p2,'DataIndex',601,'FontSize',12);


xlim([-100, 1000]);
xticks(0:100:900)
xlabel('Delay index','Position', [420 -0.08],'FontSize',12,'Interpreter','none');
ylabel('Correlation','Position', [-200 0.5], 'FontSize',12,'Interpreter','none');

set(gca,'FontSize',12)
set(gcf,'position',[600,300,700,500]);


% break y axis to two parts to show correlation clearly
breakyaxis([0.2 0.8]);

% Please click the upper part of the figure
% then run this statement in command line window
legend('Correlation between $\mathbf{X_{u1,0}}$ and $\mathbf{X_{u1,0}}$', ...
'Correlation between $\mathbf{X_{u1,0}}$ and $\mathbf{X_{u1,1}}$', ...
'Interpreter','latex','FontSize',12,'Color','none');

return;