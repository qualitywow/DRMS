function simPara = ...
    getSimParameters(scenario, Nzc, SNR, L, M, SCS, Tcp, min_cfo, max_cfo)

    simPara.fc = 30e9; % Hz, carrier  frequency, Ka band 30e9
    simPara.prachSCS = SCS; % Hz, subcarrier spacing, 60e3, 120e3

    simPara.Tseq = 1/simPara.prachSCS;

    simPara.hSat = 600e3; % 600 km
    simPara.vSat = 7.56e3; % 7.56 km/s
    simPara.Tcp = Tcp;
    simPara.min_cfo = min_cfo;
    simPara.max_cfo = max_cfo;

    simPara.cpLen = ceil((simPara.Tcp * Nzc) / (simPara.Tseq));
    simPara.numActive = L;
    simPara.numUser = M;
    simPara.SNR = SNR;
    simPara.txPower = 23; % dBm = 200 mW
end


% alpha: elevation angle
% function Tcp = getTcpByElevationAngle(hSat, minAlpha)
%     c = physconst('lightspeed'); % m/s
%     d1 = getDistance(hSat, 90); % m
%     d2 = getDistance(hSat, minAlpha);
%     Tcp = 2*(d2-d1) / c; % s
% end


% function d = getDistance(hSat, alpha)
%     R = 6371e3; % m
%     d = sqrt((R * sind(alpha))^2 + (hSat)^2 + 2*hSat*R) - (R * sind(alpha));
% end