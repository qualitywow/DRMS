function detected = ...
    detectRxSignal(approach, preambles, rxSignal, simPara, seqPara)
    detected = [];

    if (strcmp(approach, 'JSAC'))
        detected = detectJSACNew(approach, rxSignal, simPara, seqPara);
    elseif (strcmp(approach, 'ZTE'))
        detected = detectZTENew(approach, rxSignal, simPara, seqPara);
    elseif (strcmp(approach, 'HSCC'))
        detected = detectNew(approach, rxSignal, simPara, seqPara);
        % detected = detectNewN(approach, rxSignal, simPara, seqPara);
    end
end