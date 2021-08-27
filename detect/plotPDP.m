function plotPDP(Nzc, peaks, threshold, PDP, str)
    return;
    if (length(peaks) > 0)
        loc = [peaks.loc];
    else
        loc = [];
    end
    figure();
    plot(PDP, 's-'); hold on;

    % plot(ones(Nzc, 1)*prctile(PDP, 85), '-', 'LineWidth', 2); hold on;
    plot(ones(Nzc, 1)*threshold, '-', 'LineWidth', 2); hold on;
    plot(loc, PDP(loc), '^');
    hold off;
    % axis([0-30, Nzc-1+30, 0, max(PDP)+10]);
    axis([0-30, Nzc-1+30, 0, 600]);
    title(str);
    legend('PDP', 'threshold', '> threshold');
    % legend('PDP', 'prctile 80', 'threshold', '> threshold');
    % pause;
    % saveas(gcf, strcat('old-KN/',str,'.png'));
    % close();
end