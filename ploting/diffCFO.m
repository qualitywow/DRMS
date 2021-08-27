clear; clc; close all;
addpath('detect');

Nzc = 839; Ncs = 29;

CFO = [0, 0.4, 0.7, 1, 1.4, 1.7];

for i = 1:length(CFO)
    cfo = CFO(i);
    m = [];
    fprintf('CFO = %.1f\n', cfo);

    u = 45;

    % practically
    seq = zadoffChuSeq(u, Nzc);

    % s = circshift(s, 50);
    s = seq;
    out = DopplerShift(s, cfo);

    pdp = getPDP(Nzc, out, seq);

    figure();
    p = plot(0:Nzc-1, pdp, 'o-', 'LineWidth', 1.2, 'MarkerSize', 8); hold on;
    p.DataTipTemplate.DataTipRows(end).Format = '%.4f';

    a = 'Correlation under \epsilon = ';
    b = sprintf('%.1f', cfo);
    str = strcat(a, b);
    legend(str,'FontSize',16,'Color','none');
    xlim([-50, 900]);
    ylim([0, 1.05]);
    yticks([0:0.1:1]);
    xlabel('Delay index','Position', [420 -0.08],'FontSize',12,'Interpreter','none');
    ylabel('Correlation','Position', [-120 0.5], 'FontSize',12,'Interpreter','none');

    set(gca,'FontSize',12)
    set(gcf,'position',[600,300,700,500]);

    dt1 = datatip(p,'DataIndex',0+1,'FontSize',15);
    dt2 = datatip(p,'DataIndex',261+1,'FontSize',15);
    dt3 = datatip(p,'DataIndex',522+1,'FontSize',15);
    saveas(gcf, strcat('../', b, '.svg'));
end


return;