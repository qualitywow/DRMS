function outSignal = DopplerShift(signal, normCFO) %, alpha)
    % if (strcmp(simPara.scenario, 'GEO 36000'))
        % cfo = simPara.CFO;
    % else 
        % cfo = calculateCFO(simPara.vSat, simPara.hSat, simPara.fc, alpha);
    % end

    signalLen = length(signal);
    outSignal = zeros(signalLen, 1);
    
    for i = 1:signalLen
        outSignal(i) = ...
            signal(i) * exp(j*2*pi*normCFO*i / signalLen);
    end
end