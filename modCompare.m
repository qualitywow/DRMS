% clear; clc;
% OFDM modulation comparation

% nfft = 64;
% cplen = [16 32];
% nSym = 2;
% dataIn = complex(randn(nfft,nSym),randn(nfft,nSym));
% y1 = ofdmmod(dataIn,nfft,cplen);
%
% x1 = ofdmdemod(y1,nfft,cplen);
% max(x1-dataIn)



nfft = 80;
cplen = [32 0];
nullidx = [1:10 75:80].';
% nullidx = [1:6 33 64-4:64]
w = 64; h = 2;
dataIn = complex(randn(w,h),randn(w,h));
y = ofdmmod(dataIn,nfft,cplen,nullidx);
x = ofdmdemod(y,nfft,cplen,cplen,nullidx);
max(x-dataIn)

return;

Nzc = 839; K = 8; setSize = 64;
U = randi([1, Nzc-1], 3, 1);
Bcs = [29,31,37,41,43,47,53,59,61,67,71, ...
73,79,83,89,97,101,103,107,109,113,127].';
pattern = [0 1 3 6 6 3 1 0];

% setA = collectJSAC(U, Bcs, Nzc, K, setSize);
setB = collectSYMM(U, Bcs, Nzc, K, setSize, pattern);

nfft = Nzc; nifft = 1024; cplen = 180*ones(1, K);
nullidx = [1:ceil((nifft-Nzc) / 2) ceil((nifft-Nzc) / 2)+Nzc+1:nifft].';
freq = cell(setSize, 1);
for i = 1:setSize
    freq{i} = zeros(Nzc, K);
    for v = 0:K-1
        part = setB{i}(v*Nzc+1:(v+1)*Nzc);
        freq{i}(:, v+1) = fft(part, Nzc);
    end
end

y = freq{1};
yt = reshape(y, [Nzc, K]);
y_mod = ofdmmod(y, nifft, cplen, nullidx);
y_demod = ofdmdemod(y_mod, nifft, cplen, cplen, nullidx);


function preambles = collectJSAC(U, Bcs, Nzc, K, setSize)
    preambles = cell(setSize, 1);

    i = 1;
    for j = 1:length(U)
        u = U(j);
        seq = zadoffChuSeq(u, Nzc);
        for b = 1:length(Bcs)
            preamble = zeros(Nzc, 1);
            for v = 0:K-1
                preamble(v*Nzc+1:(v+1)*Nzc) = circshift(seq, -(v*Bcs(b)));
            end
            preambles{i} = preamble;
            i = i+1;

		    if (i > setSize) break; end
		end
    end
end


function preambles = collectSYMM(U, Bcs, Nzc, K, setSize, pattern)
    preambles = cell(setSize, 1);
    assert(K == length(pattern), "K is not equals to pattern size");

    i = 1;
    for j = 1:length(U)
        u = U(j);
        seq = zadoffChuSeq(u, Nzc);
        for b = 1:length(Bcs)
            preamble = zeros(Nzc, 1);
            for v = 0:K-1
                % fprintf(fid, "Cv = %d * %d\n", pattern(v+1), Bcs(b));
                preamble(v*Nzc+1:(v+1)*Nzc) = circshift(seq, -(pattern(v+1)*Bcs(b)));
            end
            preambles{i} = preamble;
            i = i+1;

            if (i > setSize) break; end
        end
    end
end


function rxSignal = passThruCh(approach, txSignals, Nzc, K, simPara, Tau, Alpha)
    rxSignal = complex(zeros(length(txSignals{1}),1), zeros(length(txSignals{1}),1));

%     ricianlChan = comm.RicianChannel('KFactor', 10);

    ricianlChan = comm.RicianChannel('KFactor', 10, ...
                                     'RandomStream','mt19937ar with seed', ...
                                     'Seed',73);

    rs = RandStream('mt19937ar','Seed',5489);

    for i = 1:length(txSignals)
        signal = circshift(txSignals{i}, Tau(i));
        if (strcmp(approach, 'DUAL-ROOT'))
            signal = phaseShiftPerSeq(signal, simPara, Nzc, K, Alpha(i));
        else
            signal = phaseShift(signal, simPara, Alpha(i));
        end

        signal = ricianlChan(signal);

        % txPower(dBm) arg(dbW), dBW = dBm - 30
%         rxSignal = rxSignal + awgn(signal, simPara.SNR, simPara.txPower-30);
        rxSignal = rxSignal + awgn(signal, simPara.SNR, simPara.txPower-30, rs);
        reset(rs);
        reset(ricianlChan); % reset channel right after use
    end

    % wgn(snr) ===
    % variance = 10^(-snr/10);
    % noise = sqrt(variance)*randn(size(x));
    % https://www.mathworks.com/matlabcentral/answers/40772-snr-in-awgn
    %{
    % tdl = nrTDLChannel('DelayProfile', 'TDL-C','FadingDistribution', 'Rician', ...
    %                'KFactor', 10.224, ...
    %                'TransmissionDirection', 'Uplink', ...
    %                'MaximumDopplerShift', simPara.CFO, ...
    %                'NumReceiveAntennas', 1);
    %                % 'DelaySpread', 100e-9, ...
    % txSignals{i} = tdl(txSignals{i});
    %}
end


function outSignal = phaseShift(signal, simPara, alpha)
    if (strcmp(simPara.scenario, 'GEO 36000'))
        cfo = simPara.CFO;
    else
        cfo = calculateCFO(simPara.vSat, simPara.hSat, simPara.fc, alpha);
    end

    signalLen = length(signal);
    outSignal = zeros(signalLen, 1);

    for i = 1:signalLen
        outSignal(i) = signal(i) * exp(j*2*pi*cfo*i / (simPara.prachSCS*signalLen));
    end
end