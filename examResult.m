function result = examResult(detected, selIdx, pTable, Nzc, simPara)
    global fid output
    sp = simPara.Tseq / Nzc; % sampling period

    result.timingMSE = 0;
    result.firstAcc = 0;
    if (length(detected) > 0)
        li = [detected.l];
        result.falseRate = length(setdiff(li, selIdx)) / length(detected);
        result.missRate = length(setdiff(selIdx, li)) / length(unique(selIdx));
        inter = intersect(li, selIdx);
        rm = []; % keep the detection that sent by UE and be detected successfully
        for lend = 1:length(detected)
            if (isempty(find(inter == li(lend))))
                rm = [rm lend];
            end
        end

        detected(rm) = [];

        % lookup the selectedPreamble table to find the actual time delay of the preamble
        sumTauI = 0;
        success = 0;
        for n = 1:length(detected)
            idx = (pTable.Preamble == detected(n).l);
            deltaTauI = abs(detected(n).tau - pTable.Tau(idx)) / Nzc;
            tauDiff = abs(detected(n).tau - pTable.Tau(idx));

            timingError = (detected(n).tau - pTable.Tau(idx)) * sp;
            fprintf(fid, "l = %2d, %f\n", detected(n).l, deltaTauI);
            sumTauI = sumTauI + deltaTauI;
            % result.timingMSE = result.timingMSE + power(deltaTauI, 2);
            % result.timingMSE = result.timingMSE + power(timingError, 2);
            result.timingMSE = result.timingMSE + power(tauDiff, 2);

            if (deltaTauI == 0) success = success + 1; end
        end
        result.avgDeltaTauI = sumTauI / length(detected);
        result.firstAcc = success / size(pTable, 1);

    else
        fprintf(fid, "No preamble signal be detected.\n");
        result.falseRate = 0.0;
        result.missRate = 1.0;
        result.avgDeltaTauI = nan;
        result.timingMSE = nan;
        result.firstAcc = 0.0;
    end


    fprintf(fid, "false alarm rate = %.4f\n", result.falseRate);
    fprintf(fid, "miss detection rate = %.4f\n", result.missRate);
    fprintf(fid, "avgDeltaTauI = %.4f\n", result.avgDeltaTauI);
    fprintf(fid, "timing MSE = %.4e\n", result.timingMSE);
    fprintf(fid, "first access rate = %.4f\n", result.firstAcc);
    outstr = sprintf('%.4f, %.4f, %.4f, %.4e, %.4f', ...
        result.falseRate, result.missRate, result.avgDeltaTauI, ...
        result.timingMSE, result.firstAcc);
    fprintf(output, '%s\n', outstr);

end

function printDetected(detected)
    for i = 1:length(detected)
        fprintf(fid, "preamble %d, tau = %d", detected(i).l, detected(i).tau);
    end
end