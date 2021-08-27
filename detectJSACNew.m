function detected = detectJSACNew(approach, rxSignal, simPara, seqPara)
    global fid;

    U = seqPara.U;
    U2 = seqPara.U;
    K = seqPara.K;
    Nzc = seqPara.Nzc;
    Bcs = seqPara.Bcs;

    detected = [];

    % get received sequence
    [rx, rxACC] = reversedProcess(rxSignal, simPara.cpLen, Nzc, K);

    % noise scale
    variance = 10^(-simPara.SNR/10);
    noise = sqrt(variance) * randn(size(rx));

    % detect root by root
    for i = 1:length(U)
        info = getFreqOffsetDist(U(i), U2(i), Nzc, K, simPara, 'KN');
        activeZC = zadoffChuSeq(U(i), Nzc);

        fprintf(fid, 'correlate rxSeq with activeZC %d\n', i);

        % pdp1 = getPDP(Nzc, rxACC, activeZC);
        % figure();
        % plot(0:Nzc-1, pdp1, '^-'); hold off;

        % split rx into K subsequences and
        % calculate the PDP of the corresponding root ZC sequence
        % then concatenate the K subsequences into pdpE
        pdpE = []; % PDP of received sequence
        pdpFA = []; % PDP of noise
        for k = 0:K-1
            temp = rx(k*Nzc+1:(k+1)*Nzc); no = noise(k*Nzc+1:(k+1)*Nzc);
            pdpE = cat(1, pdpE, getPDP(length(temp), temp, activeZC));
            pdpFA = cat(1, pdpFA, getPDP(length(no), no, activeZC));
        end
        % figure();
        % plot(0:K*Nzc-1, pdpE, 's-'); hold on;

        [threshold, peak] = ...
            setThreshold(approach, i, pdpE, pdpFA, Nzc*K, simPara, seqPara);

        % delta_1
        expect = unique(info.d);
        for p = 1:length(peak)
            if (peak(p).used == true) continue; end
            for b = 1:length(Bcs)
                % find corresponding peak group of the (root, cyclic shift)
                [idxU1, leftmostLoc1] = ...
                    peakGroup(peak, peak(p).loc, pdpE, Bcs(b), Nzc, K);
                if (K~=1 && isempty(idxU1))
                    peak(p).used = true;
                    continue;
                end

                for ii = 1:length(idxU1)
                    fprintf('%d ', peak(idxU1(ii)).loc);
                end
                fprintf('\n');
                ;

                l = convertL(Nzc, Bcs, i, Bcs(b), simPara);
                if (l == 0) continue; end
                fprintf(fid, 'detected preamble %d\n', l);

                % used
                for idx = 1:length(idxU1)
                    peak(idxU1(idx)).used = true;
                end

                % calculate fractional tau
                tauLate = mod(info.l1(end) - leftmostLoc1 + 1, Nzc);
                fprintf(fid, 'tauLate = %d\n', tauLate);


                % compensation for fractional tau and frequency offset
                rx_ref = circshift(rx, -tauLate);
                rx_ref = circshift(rx_ref, info.l1((end)));
                fprintf(fid, 'l1 = %d, left1 = %d\n', ...
                        info.l1(end), leftmostLoc1);

                C = zeros(K, 1);
                ref = circshift(activeZC, 0);
                for k = 1:K % align with 0xNcs
                    rxPart = rx_ref((k-1)*Nzc+1:k*Nzc);
                    C(k) = abs(sum((rxPart .* conj(ref)))) / Nzc;
                end


                k = find(C == max(C)) - 1; % argmax;
                fprintf(fid, 'k = %d\n', k);
                fprintf(fid, 'C(k)-----\n');
                fprintf(fid, '%f %f %f %f %f %f %f %f\n', C.');
                fprintf(fid, '\n---------');

                ;

                % bad idea but it works

                % without the second root u2,
                % JSAC can only decide the time delay by a root

                % the time delay calculation is related to ...
                % the corresponding frequency offset
                tidx = 3; % selected a specific value for JSAC

                if (K > 2)
                    if (k == 0)
                        tauLate = K*Nzc - (Nzc-tauLate);
                    else
                        if (mod(tidx,2) == 1)
                            tauLate = (k-1) * Nzc + tauLate;
                        else % (mod(tidx,2) == 0)
                            if (k == 1)
                                tauLate = k * Nzc + tauLate;
                            else
                                tauLate = (k-2) * Nzc + tauLate;
                            end
                        end
                    end
                end

                % --------------------------------------------------------------------------------
                fprintf(fid, 'tauLate: %d\n', tauLate);
                % store the detection information
                dInfo = struct('i', i, 'Ncs', Bcs(b), 'l', l, 'tau', tauLate);
                detected = [detected dInfo];
                if (length(find([detected.l] == l)) > 1)
                    ii = find([detected.l] == l, 1);
                    detected(end) = [];
                end
            end
        end
        ;
    end


end


function x = nextStep(ori, Ncs, Nzc, type)
    % which part the ori belongs to
    % the index i in equation (14)
    part = ceil(ori / Nzc);

    if (type == '+')
        if (ori+Ncs > Nzc*part)
            x = ori + Ncs;
        else
            x = ori + Ncs + Nzc;
        end
    else % type == '-'
        if (ori-Ncs < Nzc*(part-1))
            x = ori - Ncs;
        else
            x = ori - Ncs - Nzc;
        end
    end
end

function [indices, walk] = peakGroup(peak, loc0, PDP, Ncs, Nzc, K)
    indices = []; walk = loc0;

    % the possible left and right delay index for loc0
    locL = mod(nextStep(loc0, Ncs, Nzc, '-'), Nzc*K);
    locR = mod(nextStep(loc0, Ncs, Nzc, '+'), Nzc*K);

    if (isempty(find([peak.loc] == locL)) && ...
        isempty(find([peak.loc] == locR)))
        return; % can't find peak at +- Ncs => failed
    end

    n = 0; % record the number of found peak
    walk = loc0;
    while (~isempty(find([peak.loc] == mod(nextStep(walk, Ncs, Nzc, '+'), Nzc*K))) ...
           && n <= (K-2)) % keep going if there's a peak
        % record peak location
        indices = [indices find([peak.loc] == mod(nextStep(walk, Ncs, Nzc, '+'), Nzc*K))];
        n = n + 1;
        % goto next location
        walk = mod(nextStep(walk, Ncs, Nzc, '+'), Nzc*K);
    end

    walk = loc0;
    while (~isempty(find([peak.loc] == mod(nextStep(walk, Ncs, Nzc, '-'), Nzc*K))) ...
           && n <= (K-2))
        indices = [indices find([peak.loc] == mod(nextStep(walk, Ncs, Nzc, '-'), Nzc*K))];
        n = n + 1;
        walk = mod(nextStep(walk, Ncs, Nzc, '-'), Nzc*K);
    end

    if (n < (K-1)) % cannot find enough peak => failed
        indices = [];
        return;
    end

end

function [indices, walk] = peakGroupOld(peak, loc0, PDP, Ncs, Nzc, K)
    indices = []; walk = loc0;

    locL = mod(loc0 - Ncs-Nzc, Nzc*K);
    locR = mod(loc0 + Ncs+Nzc, Nzc*K);

    if (isempty(find([peak.loc] == locL)) && ...
        isempty(find([peak.loc] == locR)))
        return; % can't find peak at +- Ncs
    end

    n = 0;
    walk = loc0;
    while (~isempty(find([peak.loc] == mod(walk + Ncs + Nzc, Nzc*K))) ...
           && n <= (K-2))
        indices = [indices find([peak.loc] == mod(walk + Ncs + Nzc, Nzc*K))];
        n = n + 1;
        walk = mod(walk + Ncs + Nzc, Nzc*K);
    end

    walk = loc0;
    while (~isempty(find([peak.loc] == mod(walk - Ncs - Nzc, Nzc*K))) ...
           && n <= (K-2))
        indices = [indices find([peak.loc] == mod(walk - Ncs - Nzc, Nzc*K))];
        n = n + 1;
        walk = mod(walk - Ncs - Nzc, Nzc*K);
    end

    if (n < (K-1))
        indices = [];
        return;
    end

end


function exist = checkPeak(loc, peak)
    % global fid;
    % fprintf(fid, 'Cannot find u2 peak\n');
    if (isempty(peak))
        exist = false;
    elseif (isempty(find([peak.loc] == loc)))
        exist = false;
    else
        exist = true;
    end
end
