function detected = detectNew(approach, rxSignal, simPara, seqPara)
    global fid PPP;

    U = seqPara.U;
    U2 = seqPara.U2;
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
        activeZC2 = zadoffChuSeq(U2(i), Nzc);
        fprintf(fid, 'correlate rxSeq with activeZC %d\n', i);

        % pdp1 = getPDP(Nzc, rxACC, activeZC);
        % pdp2 = getPDP(Nzc, rxACC, activeZC2);
        % figure();
        % plot(0:Nzc-1, pdp1, '^-'); hold on;
        % plot(0:Nzc-1, pdp2, '^-'); hold off;

        % split rx into K subsequences and
        % calculate the PDP of the corresponding root ZC sequence
        % then concatenate the K subsequences into pdpE
        % HSCC use two root sequence
        pdpE = [];
        pdpE2 = []; pdpFA = [];
        for k = 0:K-1
            temp = rx(k*Nzc+1:(k+1)*Nzc); no = noise(k*Nzc+1:(k+1)*Nzc);
            pdpE = cat(1, pdpE, getPDP(length(temp), temp, activeZC));
            pdpE2 = cat(1, pdpE2, getPDP(length(temp), temp, activeZC2));
            pdpFA = cat(1, pdpFA, getPDP(length(no), no, activeZC));
        end
        PPP{i} = pdpE;
        % figure();
        % plot(0:K*Nzc-1, pdpE, 's-'); hold on;
        % plot(0:K*Nzc-1, pdpE2, 's-'); hold off;

        [threshold, peak] = ...
            setThreshold(approach, i, pdpE, pdpFA, Nzc*K, simPara, seqPara);
        [threshold, peak_u2] = ...
            setThreshold(approach, -i, pdpE2, pdpFA, Nzc*K, simPara, seqPara);

        % delta_1
        expect = unique(info.d); % d = []
        for p = 1:length(peak)
            if (peak(p).used == true) continue; end
            for b = 1:length(Bcs)
                % find corresponding peak group of the (root, cyclic shift)
                [idxU1, leftmostLoc1] = ...
                    peakGroup(peak, peak(p).loc, pdpE, Bcs(b), Nzc, K);
                if (isempty(idxU1))
                    peak(p).used = true;
                    continue;
                end

                for ii = 1:length(idxU1)
                    fprintf('%d ', peak(idxU1(ii)).loc);
                end
                fprintf(', leftmost = %d \n', leftmostLoc1);


                % find delta_2
                possibleU2_0 = [mod(leftmostLoc1 + expect, Nzc*K), ...
                               mod(leftmostLoc1 - expect, Nzc*K)]
                ;
                pass = arrayfun(@(x) checkPeak(x, peak_u2), possibleU2_0);
                pass = find(pass == true, 1);

                if (isempty(pass))
                    peak(p).used = true;
                    continue;
                end

                % the corresponding delta_2
                leftmostLoc2 = possibleU2_0(pass);
                fprintf(fid, 'u1(0): %d, u2(0): %d, Ncs = %d, pass = %d\n', ...
                       leftmostLoc1, leftmostLoc2, Bcs(b), pass);
                ;

                t1 = mod(info.l1 - leftmostLoc1, Nzc);
                t2 = mod(info.l2 - leftmostLoc2, Nzc);
                ;

                for tidx = 1:length(t1)
                    if (t1(tidx) == t2(tidx)) break; end
                end


                if (t1(tidx) ~= t2(tidx)) continue; end
                fprintf('t1(tidx) = %d, t2(tidx) = %d\n', t1(tidx), t2(tidx));

                l = convertL(Nzc, Bcs, i, Bcs(b), simPara);
                if (l == 0) continue; end
                fprintf(fid, 'detected preamble %d\n', l);


                % used
                for idx = 1:length(idxU1)
                    peak(idxU1(idx)).used = true;
                end


                % calculate fractional tau
                tauLate = mod(info.l2(tidx) - leftmostLoc2 + 1, Nzc);
                fprintf(fid, 'tauLate f: %d\n', tauLate);

                % compensation for fractional tau and frequency offset
                rx_ref = circshift(rx, -tauLate);
                rx_ref = circshift(rx_ref, info.l2((tidx)));
                ;


                C = zeros(K, 1);
                ref = circshift(activeZC2, 0);
                for k = 1:K % align with u2
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

                % the time delay calculation is related to ...
                % the corresponding frequency offset

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
                fprintf(fid, 'tauLate i: %d\n', tauLate);
                % store the detection information
                dInfo = struct('i', i, 'Ncs', Bcs(b), 'l', l, 'tau', tauLate, 'maxC', max(C));
                detected = [detected dInfo];
                if (length(find([detected.l] == l)) > 1)
                    ii = find([detected.l] == l, 1);
                    if (detected(ii).maxC > max(C))
                        detected(end) = [];
                    else
                        detected(ii) = [];
                    end
                end
            end
        end
        % plotPDP(Nzc*K, peaks, threshold, pdpE, 'pdpE');
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

    locL = mod(nextStep(loc0, Ncs, Nzc, '-'), Nzc*K);
    locR = mod(nextStep(loc0, Ncs, Nzc, '+'), Nzc*K);

    if (isempty(find([peak.loc] == locL)) && ...
        isempty(find([peak.loc] == locR)))
        return; % can't find peak at +- Ncs
    end

    n = 0;
    walk = loc0;
    while (~isempty(find([peak.loc] == mod(nextStep(walk, Ncs, Nzc, '+'), Nzc*K))) ...
           && n <= (K-2))
        indices = [indices find([peak.loc] == mod(nextStep(walk, Ncs, Nzc, '+'), Nzc*K))];
        n = n + 1;
        walk = mod(nextStep(walk, Ncs, Nzc, '+'), Nzc*K);
    end

    walk = loc0;
    while (~isempty(find([peak.loc] == mod(nextStep(walk, Ncs, Nzc, '-'), Nzc*K))) ...
           && n <= (K-2))
        indices = [indices find([peak.loc] == mod(nextStep(walk, Ncs, Nzc, '-'), Nzc*K))];
        n = n + 1;
        walk = mod(nextStep(walk, Ncs, Nzc, '-'), Nzc*K);
    end

    if (n < (K-1))
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
