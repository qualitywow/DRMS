function info = getFreqOffsetDist(u1, u2, Nzc, K, simPara, type)
    info = struct();
    if (strcmp(type, 'N')) % ZTE derived
        info = foffsetN(u1, u2, Nzc, K, simPara);
    elseif (strcmp(type, 'KN')) % HSCC derived
        info = foffsetKN(u1, u2, Nzc, K, simPara);
    end
end

function info = foffsetN(u1, u2, Nzc, K, simPara)
    % apply extexded euclidean algorithm to calculate the shift offset of u1, u2
    % need the background knowlegde of ZTE paper
    global fid;
    info = struct();
    [x1, y1, gcd1] = extended_euclid(u1, Nzc);
    [x2, y2, gcd2] = extended_euclid(u2, Nzc);
    if (gcd1 ~= 1 || gcd2 ~= 1)
        fprintf(fid, "check for freq distance\n");
    end

    info.s1 = Nzc - x1;
    info.s2 = Nzc - x2;
    roundEpsilon = round(epsilon);
    % roundEpsilon = round(epsilon/K) then the answer is exactly what I derived

    info.l1 = mod(info.s1*roundEpsilon, Nzc);
    info.l2 = mod(info.s2*roundEpsilon, Nzc);
    info.d = mod(roundEpsilon * (info.s1-info.s2), Nzc);

    fprintf(fid, "========================\n");
    fprintf(fid, "u1 = %d, u2 = %d\n", u1, u2);
    fprintf(fid, "s1 = %d, s2 = %d\n", info.s1, info.s2);
    fprintf(fid, "l1 = %d, l2= %d, dist = %d\n", info.l1, info.l2, info.d);
end

function [x,y,d]=extended_euclid(a,b)
    if b==0
        x=1;
        y=0;
        d=a;
        return;
    end
    [x1,y1,d1]=extended_euclid(b,mod(a,b));

    x=y1;
    y=x1-floor(a/b)*y1;
    d=d1;

end


function info = foffsetKN(u1, u2, Nzc, K, simPara)
    global fid;
    info = struct();
    if (simPara.prachSCS == 60e3 || simPara.prachSCS == 120e3)
        ar = simPara.min_cfo:simPara.max_cfo;
        range = sort(unique(round(ar / K)), 'descend');
        roundEpsilon = [];
        for r = 1:length(range)*2
            if (mod(r,2) == 0)
                signed = -1;
            else
                signed = 1;
            end
            roundEpsilon = [roundEpsilon range(ceil(r/2))*signed];
        end
    % elseif (simPara.prachSCS == 5e3)
    %     % assume round epsilon = 10
    %     roundEpsilon = [round(10/K), round(10/K)* (-1)];
    % elseif (simPara.prachSCS == 1.25e3)
    %     roundEpsilon = [round(40/K), round(41/K)];
    end

    % given specific CFO
    % roundEpsilon = round(simPara.CFO/simPara.prachSCS/K);

    x1 = SolveCongruence(u1, roundEpsilon, Nzc);
    x2 = SolveCongruence(u2, roundEpsilon, Nzc);

    info.roundEpsilon = roundEpsilon;
    info.l1 = Nzc - x1;
    info.l2 = Nzc - x2;

    info.d = abs(info.l1 - info.l2);
    % info.d = mod(info.l1 - info.l2, Nzc);

    fprintf(fid, "========================\n");
    fprintf(fid, "u1 = %d, u2 = %d\n", u1, u2);
    for i = 1:length(roundEpsilon)
        fprintf(fid, "l1 = %d, l2= %d, dist = %d\n", info.l1(i), info.l2(i), info.d(i));
    end
end