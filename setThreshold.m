function [threshold, peaks] = setThreshold(approach, zcIdx, PDP, PDP_FA, Nzc, simPara, seqPara)
    % [maxVal, minVal] = maxPeakNum(approach, simPara, seqPara);
    if (strcmp(approach, 'JSAC'))
        threshold = setJSAC(approach, zcIdx, PDP, PDP_FA, Nzc, simPara, seqPara);
    elseif (strcmp(approach, 'ZTE'))
        threshold = setZTE(approach, zcIdx, PDP, PDP_FA, Nzc, simPara, seqPara);
    elseif (strcmp(approach, 'HSCC'))
        threshold = setHSCC(approach, zcIdx, PDP, PDP_FA, Nzc, simPara, seqPara);
        % accumulating detection
        % threshold = setHSCCN(approach, zcIdx, PDP, PDP_FA, Nzc, simPara, seqPara);
    end

    disp(mean(PDP) / mean(PDP_FA));

    loc = find(PDP > threshold);
    if (1) % length(loc) > 0
        if (length(loc) > 0)
            for p = 1:length(loc)
                peaks(p).loc = loc(p);
                peaks(p).used = false;
            end
        else
            peaks = [];
        end

        str = sprintf('%s- PDP of rxSignal and activeZC{%d}', approach, zcIdx);
        plotPDP(Nzc, peaks, threshold, PDP, str);
        % fprintf(fid, "numOfPeak: %d\n", length(peaks));
        % fprintf(fid, '%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n', ...
                % mean(PDP), max(PDP), min(PDP), ...
                % prctile(PDP, 25), prctile(PDP, 50), prctile(PDP, 75));
        ;
    else
        peaks = [];
        % fprintf(fid, "numOfPeak: %d\n", 0);
    end
end

% the eta value is helping to set a threshold that ...
% >>> make the FALSE ALARM RATE = 0.0. <<<
% please adjust the value if the simulation parameter is changed,
% but don't delete current values, or you need to find it through the simulations
function threshold = setJSAC(approach, zcIdx, PDP, PDP_FA, Nzc, simPara, seqPara)
    global fid
    etaK4 = [4, 10, 20];
    etaK8 = [-4, -2.5, 0.5];
    if (seqPara.K == 4)
        etaK = etaK4;
    elseif (seqPara.K == 8)
        etaK = etaK8;
    end

    if (mean(PDP) / mean(PDP_FA) < 2000)
        eta = etaK(1);
    elseif (mean(PDP) / mean(PDP_FA) < 4000)
        eta = etaK(2);
    elseif (mean(PDP) / mean(PDP_FA) < 10000)
        eta = etaK(3);
    else
        eta = 10;
    end
    threshold = (max(PDP_FA) / mean(PDP_FA)+eta) * mean(PDP_FA) * Nzc;

end

function threshold = setZTE(approach, zcIdx, PDP, PDP_FA, Nzc, simPara, seqPara)
    global fid
    etaK4 = [13, 23, 23];
    etaK8 = [-2.5, 0.5, 4.5];
    if (seqPara.K == 4)
        etaK = etaK4;
    elseif (seqPara.K == 8)
        etaK = etaK8;
    end

    if (mean(PDP) / mean(PDP_FA) < 2000)
        eta = etaK(1);
    elseif (mean(PDP) / mean(PDP_FA) < 4000)
        eta = etaK(2);
    elseif (mean(PDP) / mean(PDP_FA) < 10000)
        eta = etaK(3);
    else
        eta = 10;
    end
    threshold = (max(PDP_FA) / mean(PDP_FA)+eta) * mean(PDP_FA) * Nzc;

end

function threshold = setHSCC(approach, zcIdx, PDP, PDP_FA, Nzc, simPara, seqPara)
    global fid
    etaK4 = [1, 8, 21];
    etaK8 = [-4.5, -2, 0];
    if (seqPara.K == 4)
        etaK = etaK4;
    elseif (seqPara.K == 8)
        etaK = etaK8;
    end

    if (mean(PDP) / mean(PDP_FA) < 2000)
        eta = etaK(1);
    elseif (mean(PDP) / mean(PDP_FA) < 4000)
        eta = etaK(2);
    elseif (mean(PDP) / mean(PDP_FA) < 10000)
        eta = etaK(3);
    else
        eta = 10;
    end
    threshold = (max(PDP_FA) / mean(PDP_FA)+eta) * mean(PDP_FA) * Nzc;
end


function threshold = setHSCCN(approach, zcIdx, PDP, PDP_FA, Nzc, simPara, seqPara)
    % NEED TUNING FOR SPECIFIC PDP RANGE
    global fid
    etaK4 = [80, 90, 90];
    etaK8 = [-4.5, -2, 0];
    if (seqPara.K == 4)
        etaK = etaK4;
    elseif (seqPara.K == 8)
        etaK = etaK8;
    end

    if (mean(PDP) / mean(PDP_FA) < 4000)
        eta = etaK(1);
    elseif (mean(PDP) / mean(PDP_FA) < 6000)
        eta = etaK(2);
    elseif (mean(PDP) / mean(PDP_FA) < 10000)
        eta = etaK(3);
    else
        eta = 10;
    end
    threshold = (max(PDP_FA) / mean(PDP_FA)+eta) * mean(PDP_FA) * Nzc;
 end




function [maxVal, minVal] = maxPeakNum(approach, simPara, seqPara)
    maxVal = -1;
    avgPerRoot = ceil(simPara.numUser / length(seqPara.U));
    if (strcmp(approach, 'JSAC'))
        seqPerRoot = ceil(simPara.numActive/length(seqPara.U));
        if (simPara.numUser <= seqPerRoot) % 2, 4, 8, 16
            maxVal = simPara.numUser * 8;
        else % 32, 64
            maxVal = seqPerRoot * 8;
        end
        maxVal = avgPerRoot * 8;
        minVal = 8;
    elseif (strcmp(approach(1:5), 'SYMM-'))
        maxVal = avgPerRoot * 4;
        minVal = 4;
    elseif (strcmp(approach(1:4), 'DUAL'))
        seqPerRoot = ceil(simPara.numActive/length(seqPara.U));
        if (simPara.numUser <= seqPerRoot) % 2, 4, 8, 16
            maxVal = simPara.numUser * 4;
        else % 32, 64
            maxVal = seqPerRoot * 4;
        end
        maxVal = avgPerRoot * 4;
        minVal = 4;
    end
end