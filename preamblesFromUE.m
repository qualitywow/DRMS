function [rxSignal, selectedIndices, pTable] = ...
    preamblesFromUE(approach, preambles, simPara, seqPara)

    U = seqPara.U;
    K = seqPara.K;
    Nzc = seqPara.Nzc;
    Bcs = seqPara.Bcs;

    % UE select preamble index to send (arranged)
    selectedIndices = samplePreamble(simPara.numActive, simPara.numUser, U);
    selectedPreambles = preambles(selectedIndices);

    ueInfo = getUE(simPara);

    % preambles -> signals
    txSignals = collectSignals(selectedPreambles, simPara, Nzc, K);
    rxSignal = passThruCh(approach, txSignals, simPara, Nzc, K, ueInfo);

    pTable = showSelectedPreambles(approach, U, Bcs, selectedIndices, ...
                                   ueInfo, simPara);
end

function rxSignal = passThruCh(approach, txSignals, simPara, Nzc, K, ueInfo)
    txLen = length(txSignals{1}); % [CP SEQ]
    len = txLen + simPara.cpLen; %  [CP SEQ GT]
    rxSignal = complex(zeros(len,1), zeros(len,1));

    Epsilon = [ueInfo.eps];
    Tau = [ueInfo.tau];
    SNR = simPara.SNR;

    tdl = nrTDLChannel('DelayProfile', 'Custom', ...
                       'FadingDistribution', 'Rician', ...
                       'KFactorFirstTap', 10.224, ...
                       'TransmissionDirection', 'Uplink', ...
                       'NumReceiveAntennas', 1, ...
                       'RandomStream', "mt19937ar with seed", ...
                       'Seed', 84); % NTN-TDL-C channel

    for i = 1:length(txSignals)
        signal = txSignals{i}; % [CP SEQ]

        signal = dopplerShift(signal, Epsilon(i));
        % [delay CP SEQ GT]
        signal = cat(1, zeros(Tau(i), 1), signal, ...
                     zeros(len-(Tau(i)+txLen), 1));

        signal = tdl(signal);

        variance = 10^(-SNR/10);
        noise = sqrt(variance) * randn(size(signal));
        signal = signal * simPara.txPower + noise;

        rxSignal = rxSignal + signal;
    end
end


function info = getUE(simPara)
    % rng(15); % randi random seed reset
    neg = randi([1, simPara.numUser], floor(simPara.numUser/2), 1);
    signed = ones(simPara.numUser, 1); signed(neg) = -1;

    % rng(15); % randi random seed reset
    if (simPara.prachSCS == 60e3 || simPara.prachSCS == 120e3)
        info.eps = [randi([round(simPara.min_cfo)*10 round(simPara.max_cfo)*10+4], ...
                           simPara.numUser, 1) / 10 .* signed];

        % info.eps = ones(simPara.numUser, 1) * 2.5;
        % info.eps = [1].';
    % elseif (simPara.prachSCS == 1.25e3)
    %     info.eps = randi([405, 414], simPara.numUser, 1) / 10;
    end

    % rng(15); % randi random seed reset
    info.tau = randi([0, simPara.cpLen], simPara.numUser, 1);
    % info.tau = ones(simPara.numUser, 1)*279;
end


function outSignal = dopplerShift(signal, cfo)
    signalLen = length(signal);
    outSignal = zeros(signalLen, 1);

    for i = 1:signalLen
        outSignal(i) = ...
            signal(i) * exp(j*2*pi*cfo*i / signalLen);
    end
end



function txSignals = collectSignals(selectedPreambles, simPara, Nzc, K)
    txSignals = cell(length(selectedPreambles), 1);

    for i = 1:length(selectedPreambles)
        x = selectedPreambles{i}; s = [];
        for k = 1:K
            temp = x((k-1)*Nzc+1:k*Nzc);
            temp = fft(temp); % freq-domain
            temp = ifft(temp); % time-domain
            s = cat(1, s, temp);
        end

        % tx signal = [CP SEQ]
        txSignal = cat(1, s(length(s)-simPara.cpLen+1:end), s);
        % txSignal = txSignal * simPara.txPower;
        txSignals{i} = txSignal;
    end
end

function selectedIndices = samplePreamble(numActive, numUser, U)
    rng();
    selectedIndices = randperm(numActive, numUser);
    selectedIndices = reshape(selectedIndices, [numUser, 1]);
    % selectedIndices = [];
    % m = 0; n = 1; i = 0;
    % while (i < numUser)
    %     i = i+1;
    %     selectedIndices = [selectedIndices 22*m+n];

    %     m = m+1;
    %     if (m == length(U))
    %         m = 0;
    %         n = n+1;
    %     end
    % end

    % if (~isempty(find(selectedIndices == 65)))
    %     selectedIndices(find(selectedIndices == 65)) = 44;
    % end
    % selectedIndices = selectedIndices.';
end