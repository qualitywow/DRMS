function preambles = collectPreambles(approach, setSize, seqPara)
    K = seqPara.K;
    Nzc = seqPara.Nzc;
    Bcs = seqPara.Bcs;

    if (strcmp(approach, 'JSAC'))
        preambles = collectJSAC(Bcs, Nzc, K, setSize, seqPara);
    elseif (strcmp(approach, 'ZTE'))
        preambles = collectZTE(Bcs, Nzc, K, setSize, seqPara);
    elseif (strcmp(approach, 'HSCC'))
        preambles = collectHSCC(Bcs, Nzc, K, setSize, seqPara);
    end
end

function preambles = collectJSAC(Bcs, Nzc, K, setSize, seqPara)
    %|-----------------------------------------------|
    %|  0  |  1  |  2  |  3  |  4  |  5  |  6  |  7  |
    %|-----------------------------------------------|
    %|                      u1                       |
    %|-----------------------------------------------|

    preambles = cell(setSize, 1);
    U = seqPara.U;

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


function preambles = collectZTE(Bcs, Nzc, K, setSize, seqPara)
    %|-----------------------------------------------|
    %|  0  |  1  |  2  |  3  |  0  |  1  |  2  |  3  |
    %|-----------------------------------------------|
    %|           u1          |          u2           |
    %|-----------------------------------------------|
    % 改成u2(1234), 否則K = 2會無法分出Ncs
    preambles = cell(setSize, 1);
    U = seqPara.U;
    U2 = seqPara.U2;

    i = 1;
    for j = 1:length(U)
        seq = zadoffChuSeq(U(j), Nzc);
        seq2 = zadoffChuSeq(U2(j), Nzc);
        for b = 1:length(Bcs)
            preamble = zeros(Nzc, 1);
            for v = 0:K-1
                if (v < K/2)
                    preamble(v*Nzc+1:(v+1)*Nzc) = circshift(seq, -(v*Bcs(b)));
                else
                    k = v - K/2;
                    preamble(v*Nzc+1:(v+1)*Nzc) = circshift(seq2, -(k*Bcs(b)));
                end
            end
            preambles{i} = preamble;
            i = i+1;

            if (i > setSize) break; end
        end
    end
end


function preambles = collectHSCC(Bcs, Nzc, K, setSize, seqPara)
    %|-----------------------------------------------|
    %|  0  |  1  |  2  |  3  |  4  |  5  |  6  |  7  |
    %|-----------------------------------------------|
    %|u1+u2|                u1                       |
    %|-----------------------------------------------|

    preambles = cell(setSize, 1);
    U = seqPara.U;
    U2 = seqPara.U2;

    i = 1;
    for j = 1:length(U)
        seq = zadoffChuSeq(U(j), Nzc);
        seq2 = zadoffChuSeq(U2(j), Nzc);
        for b = 1:length(Bcs)
            preamble = zeros(Nzc, 1);
            for v = 0:K-1
                if (v == 0)
                    preamble(v*Nzc+1:(v+1)*Nzc) = seq + seq2;
                else
                    preamble(v*Nzc+1:(v+1)*Nzc) = circshift(seq, -(v*Bcs(b)));
                end
            end
            preambles{i} = preamble;
            i = i+1;

            if (i > setSize) break; end
        end
    end
end