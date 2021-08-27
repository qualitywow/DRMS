function preambleTable = ...
    showSelectedPreambles(approach, U, Bcs, selectedIndices, ueInfo, simPara)
    global fid

    preambleTable = [];

    Epsilon = ueInfo.eps;
    Tau = ueInfo.tau;

    User = [1:simPara.numUser].';
    Preamble = selectedIndices;

    RootIndex = zeros(simPara.numUser, 1);
    CyclicShift = zeros(simPara.numUser, 1);

    for p = 1:length(selectedIndices)
        for i = 1:length(U)
            if (selectedIndices(p) <= (length(Bcs)*i))
                RootIndex(p) = U(i);
                break;
            end
        end

        CyclicShift(p) = mod(selectedIndices(p), length(Bcs));
        if (CyclicShift(p) == 0)
            CyclicShift(p) = Bcs(end);
        else
            CyclicShift(p) = Bcs(CyclicShift(p));
        end
    end

    preambleTable = table(User, Preamble, RootIndex, CyclicShift, Tau, Epsilon);

    fprintf(fid, 'UE\t#l\t  u\t Ncs\t Tau\t CFO\n');
    t = preambleTable.Variables;
    for row = 1:size(t, 1)
        fprintf(fid, "%2d\t%2d\t%3d\t %3d\t%4d\t%4.1f\n", t(row,:));
    end
end