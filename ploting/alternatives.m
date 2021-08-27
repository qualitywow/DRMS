clear; clc; close all;
addpath('detect');

Nzc = 139; Ncs = 29; u = 45; K = 8;

seq = zadoffChuSeq(u, Nzc);
sp = [0 1 2 3 4 5 6 7]; s = [];
for v = 0:K-1
    if (v == 0)
        s  = cat(1, s, seq); % +seq2
    else
        s  = cat(1, s, circshift(seq, -sp(v+1)*Ncs));
    end
end


out1 = DopplerShift(circshift(s, 100), 0.4);
out2 = DopplerShift(circshift(s, 200), 1.2);

pdp = getPDP(Nzc, out1+out2, seq);

figure();
p = plot(0:Nzc-1, pdp, 'o-', 'LineWidth', 1.2, 'MarkerSize', 8); hold on;
p.DataTipTemplate.DataTipRows(end).Format = '%.4f';

% a = 'Correlation under \epsilon = ';
% b = sprintf('%.1f', cfo);
% str = strcat(a, b);
% legend(str,'FontSize',16,'Color','none');
xlim([-50, 900]);
ylim([0, 1.05]);
yticks([0:0.1:1]);
xlabel('Delay index','Position', [420 -0.08],'FontSize',12,'Interpreter','none');
ylabel('Correlation','Position', [-120 0.5], 'FontSize',12,'Interpreter','none');

set(gca,'FontSize',12)
set(gcf,'position',[600,300,700,500]);

% dt1 = datatip(p,'DataIndex',0+1,'FontSize',15);
% dt2 = datatip(p,'DataIndex',261+1,'FontSize',15);
% dt3 = datatip(p,'DataIndex',522+1,'FontSize',15);
% saveas(gcf, strcat('../', b, '.svg'));

return;