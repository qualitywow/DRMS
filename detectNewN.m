% accumulating detect
% NOT USING IN THE THESIS
function detected = detectNewN(approach, rxSignal, simPara, seqPara)
    global fid PPP;
    N = 2;
    U = seqPara.U;
    U2 = seqPara.U2;
    K = seqPara.K;
    Nzc = seqPara.Nzc;
    Bcs = seqPara.Bcs;

    detected = [];


    [rx, rxACC] = reversedProcess(rxSignal, simPara.cpLen, Nzc, K);
    [rx, rxList] = reversedProcessN(rxSignal, simPara.cpLen, Nzc, K, 1);

    variance = 10^(-simPara.SNR/10);
    noise = sqrt(variance) * randn(Nzc, 1);

    for i = 1:length(U)
        info = getFreqOffsetDist(U(i), U2(i), Nzc, K, simPara, 'KN');
        activeZC = zadoffChuSeq(U(i), Nzc);
        activeZC2 = zadoffChuSeq(U2(i), Nzc);
        fprintf(fid, 'correlate rxSeq with activeZC %d\n', i);

        pdpStruct = [];
        for rxIdx = 1:length(rxList)
            pdp1 = getPDP(Nzc, rxList{rxIdx}, activeZC);
            pdp2 = getPDP(Nzc, rxList{rxIdx}, activeZC2);
            temp.pdp1 = pdp1;
            temp.pdp2 = pdp2;
            pdpStruct = [pdpStruct temp];
            % figure();
            % plot(0:Nzc-1, pdp1, '^-'); hold on;
            % plot(0:Nzc-1, pdp2, '^-'); hold off;
        end
        ;

        pdpFA = getPDP(Nzc, noise, activeZC);

        tempU1 = [];
        tempU2 = [];
        for rxIdx = 1:length(pdpStruct)

            [threshold, peak] = ...
                setThreshold(approach, i, pdpStruct(rxIdx).pdp1, pdpFA, Nzc, simPara, seqPara);
            [threshold, peak_u2] = ...
                setThreshold(approach, -i, pdpStruct(rxIdx).pdp2, pdpFA, Nzc, simPara, seqPara);

            if (~isempty(peak)) tempU1 = union(tempU1, [peak.loc]); end
            if (~isempty(peak_u2)) tempU2 = union(tempU2, [peak_u2.loc]); end
        end


        peak = []; peak_u2 = [];

        for idxTemp = 1:length(tempU1)
            peak = [peak struct('loc',tempU1(idxTemp),'used',false)];
        end
        for idxTemp = 1:length(tempU2)
            peak_u2 = [peak struct('loc',tempU2(idxTemp),'used',false)];
        end
        ;

        expect = unique(info.d); % d = []
        for p = 1:length(peak)
            if (peak(p).used == true) continue; end
            for b = 1:length(Bcs)
                % if (peak(p).used == true) break; end

%                 [idxU1, leftmostLoc1] = ...
%                     peakGroup(peak, peak(p).loc, pdpE, Bcs(b), Nzc, K);

                [idxU1, leftmostLoc1] = peakGroupOld(peak, peak(p).loc, pdp1, Bcs(b), Nzc, K);

                if (isempty(idxU1))
                    peak(p).used = true;
                    continue;
                end

                for ii = 1:length(idxU1)
                    fprintf('%d ', peak(idxU1(ii)).loc);
                end
                fprintf(', leftmost = %d \n', leftmostLoc1);

                % possibleU2_0 = [mod(peak(p).loc + expect, Nzc), ...
                %               mod(peak(p).loc - expect, Nzc)];

                possibleU2_0 = [mod(leftmostLoc1 + expect, Nzc), ...
                               mod(leftmostLoc1 - expect, Nzc)]
                ;
                pass = arrayfun(@(x) checkPeak(x, peak_u2), possibleU2_0);
                pass = find(pass == true, 1);

                if (isempty(pass))
                    peak(p).used = true;
                    continue;
                end

                leftmostLoc2 = possibleU2_0(pass);
                fprintf(fid, 'u1(0): %d, u2(0): %d, Ncs = %d, pass = %d\n', ...
                       leftmostLoc1, leftmostLoc2, Bcs(b), pass);
                ;
                % if (tauLate == Nzc) tauLate = 0; end
                t1 = mod(info.l1 - leftmostLoc1, Nzc);
                t2 = mod(info.l2 - leftmostLoc2, Nzc);
                ;

                for tidx = 1:length(t1)
                    if (t1(tidx) == t2(tidx)) break; end
                end

                % if (t1 ~= t2) continue; end
                if (t1(tidx) ~= t2(tidx)) continue; end
                fprintf('t1(tidx) = %d, t2(tidx) = %d\n', t1(tidx), t2(tidx));

                l = convertL(Nzc, Bcs, i, Bcs(b), simPara);
                if (l == 0) continue; end
                fprintf(fid, 'detected preamble %d\n', l);


                % used
                for idx = 1:length(idxU1)
                    peak(idxU1(idx)).used = true;
                end


                % calculate fractional time delay
                tauLate = mod(info.l2(tidx) - leftmostLoc2 + 1, Nzc);
                fprintf(fid, 'tauLate f: %d\n', tauLate);

%                 r = floor((info.l2(tidx) * info.roundEpsilon(tidx)) / Nzc);
%                 q = floor((leftmostLoc2-1) / Nzc);
%                 fprintf(fid, 'l2: %d, Y[%d], leftmostLoc2: %d, Y[%d]\n', ...
%                         info.l2(tidx), r, leftmostLoc2, q) ;

                %%{
                % tauLate
                rx_ref = circshift(rx, -tauLate);
                rx_ref = circshift(rx_ref, info.l2((tidx)));
                ;


                C = zeros(K, 1);
                ref = circshift(activeZC2, 0);
                for k = 1:K % 對齊 1xNcs
                    rxPart = rx_ref((k-1)*Nzc+1:k*Nzc);
                    C(k) = abs(sum((rxPart .* conj(ref)))) / Nzc;
                end


                k = find(C == max(C)) - 1; % argmax;
                fprintf(fid, 'k = %d\n', k);
                fprintf(fid, 'C(k)-----\n');
                fprintf(fid, '%f %f %f %f %f %f %f %f\n', C.');
                fprintf(fid, '\n---------');

                ;

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

                %%}

                % tauLate = tauLate + Nzc * (q);
                % --------------------------------------------------------------------------------
                fprintf(fid, 'tauLate i: %d\n', tauLate);
                ;
                %
                dInfo = struct('i', i, 'Ncs', Bcs(b), 'l', l, 'tau', tauLate, 'maxC', max(C));
                detected = [detected dInfo]; % 直接先放
                if (length(find([detected.l] == l)) > 1) %可能會有兩個
                    ii = find([detected.l] == l, 1);
                    if (detected(ii).maxC > max(C)) % 留max(C)大的
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
    % 目前ori屬於哪一段
    part = ceil(ori / Nzc);
    % 往左往右Ncs就超過目前屬於的這一段，那就不用加減Nzc讓他跳到下一段

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

    locL = mod(loc0 - Ncs, Nzc);
    locR = mod(loc0 + Ncs, Nzc);

    if (isempty(find([peak.loc] == locL)) && ...
        isempty(find([peak.loc] == locR)))
        return; % can't find peak at +- Ncs
    end

    n = 0;
    walk = loc0;
    while (~isempty(find([peak.loc] == mod(walk + Ncs, Nzc))) ...
           && n <= (K-2))
        indices = [indices find([peak.loc] == mod(walk + Ncs, Nzc))];
        n = n + 1;
        walk = mod(walk + Ncs, Nzc);
    end

    walk = loc0;
    while (~isempty(find([peak.loc] == mod(walk - Ncs, Nzc))) ...
           && n <= (K-2))
        indices = [indices find([peak.loc] == mod(walk - Ncs, Nzc))];
        n = n + 1;
        walk = mod(walk - Ncs, Nzc);
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


function [rx, rxList] = reversedProcessN(rxSignal, cpLen, Nzc, K, alpha)
    rx = zeros(K*Nzc, 1);
    rxList = cell(ceil(K/alpha), 1);
    rxACC = zeros(Nzc, 1);
    rxSignal = rxSignal(cpLen+1:end);

    ii = 1;
    ij = 0;
    count = 1;
    while ii <= K
        temp = rxSignal((ii-1)*Nzc+1:ii*Nzc);
        rxFFT = fft(temp, Nzc);
        rxSeq = ifft(rxFFT, Nzc);
        rx((ii-1)*Nzc+1:ii*Nzc) = rxSeq;

        rxACC = rxACC + rxSeq;
        ij = ij + 1;

        if (ij == alpha || ii == K)
            rxList{count} = rxACC;
            rxACC = zeros(Nzc, 1);
            ij = 0;
            count = count + 1;
        end
        ii = ii + 1;
    end
end