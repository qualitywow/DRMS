function plotResult(appName, snr, numUE, result)
    ;
    % result that same numUE under different SNR
    names = ["fa", "md", "te", "ac"];
    marker = ['s', '^', '*'];
    xticks = snr;
    for i = 1:length(names)
        plotSth(names(i), result.(names(i)), appName, marker, snr, numUE);
    end
end


function plotSth(name, sth, appName, marker, snr, numUE)
    figure();
    ax = axes;
    FaceOrder = ax.ColorOrder;
    markSize = [8, 6, 7];
    for i = 1:length(appName)
        y = cell2mat(sth(:,i));
        plot(snr, y, strcat(marker(i), "-"), ...
            'MarkerSize', markSize(i), ...
            'LineWidth', 1);
            %  'MarkerFaceColor', FaceOrder(i, :),
        hold on;
    end
    hold off;
    grid on;
    legend(appName, 'Location', 'best');

    if (strcmp(name, "fa"))
        ylim([0, 1]);
        str = sprintf("False alarm rate, numUE = %d", numUE);
    elseif (strcmp(name, "md"))
        ylim([0, 1]);
        str = sprintf("Miss detection rate, numUE = %d", numUE);
    elseif (strcmp(name, "te"))
        str = sprintf("Timing MSE, numUE = %d", numUE);
    elseif (strcmp(name, "ac"))
        ylim([0, 1]);
        str = sprintf("First access rate, numUE = %d", numUE);
    end
    title(str);
%     saveas(gcf, strcat('hardcode_threshold/', str, '.png'));
%     saveas(gcf, strcat('A/', str, '.png'));
end