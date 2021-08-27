function seqPara = getSeqParameters(scenario, approach, K)
    seqPara.Nzc = 139;

    seqPara.Bcs = [29,31,37,41,43,47,53,59,61,67,71, ...
                   73,79,83,89,97,101,103,107,109,113,127].';

    seqPara.U = [48, 55, 73]; disp(seqPara.U);
    if (strcmp(approach, 'ZTE') || strcmp(approach, 'HSCC'))
        seqPara.U2 = 50 + seqPara.U; disp(seqPara.U2);
    end

    seqPara.K = K;
end


