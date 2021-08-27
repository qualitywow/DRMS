clear; clc; close all;

addpath('detect'); % include matlab code in './detect'

% output log file
global fid;
global output;
% fid = fopen('log.txt', 'w');
fid = 1; % 1 for console output


numL = 64; % available preamble
scenario = 'LEO 600';

% approach can be one of {'JSAC', 'ZTE', 'HSCC'}
appName{1} = 'JSAC';
appName{2} = 'HSCC';
appName{3} = 'ZTE';
% appName{1} = 'HSCC';
% appName{1} = 'ZTE';

% load predefined parameters
load 'tcp_table.mat';
load 'cfo_table.mat';

% set 
r = 25; scs = 60e3; beam = 2; K = 4;

cp_idx = find((tcp_table.Radius == r) & ...
              (tcp_table.SCS == scs) & ...
              (tcp_table.Beam == beam) & ...
              (tcp_table.K == K));
Tcp = tcp_table.Tcp(cp_idx);

cfo_idx = find((cfo_table.Radius == r) & ...
               (cfo_table.SCS == scs) & ...
               (cfo_table.Beam == beam));
min_cfo = cfo_table.Min(cfo_idx);
max_cfo = cfo_table.Max(cfo_idx);

snr = [-8]; % -20:2:
ue = [1]; % 2, 4, 8, 16, 32

for n = 1:length(ue)
    UE = ue(n);
    for m = 1:length(snr)
        SNR = snr(m);
        for k = 1:length(appName)
            % output filepath
            folder = '../rawData';
            
            % namestyle
            % 1: Variable, 2: Approach, others
            filename = sprintf('snr%d_%s_r%d_%dkHz_b%d_K%d_ue%02d.csv', ...
                       SNR, appName{k}, r, round(scs/1e3, 0), beam, K, ue(n));

            % filename = sprintf('b%d_%s_snr%d_r%d_%dkHz_K%d_ue%02d.csv', ...
            % beam, appName{k}, SNR, r, round(scs/1e3, 0), K, ue(n));

            % filename = sprintf('ue_%02d_%s_snr%d_b%d_r%d_%dkHz_K%d.csv', ...
            % ue(n), appName{k}, SNR, beam, r, round(scs/1e3, 0), K);
            
            %filename = sprintf('eps%d_ue_%02d_%s_snr%d_r%d_%dkHz_K%d.csv', ...
            %2.5, ue(n), appName{k}, SNR, r, round(scs/1e3, 0), K);
            
            output = fopen(strcat(folder, "/", filename), 'a+');
            % output = 1;

            for t = 1:1
                approach = appName{k};
                fprintf(fid, '==================================\n');
                fprintf(fid, 'UE = %d, SNR = %d, design = %s\n', UE, SNR, approach);
                fprintf(fid, '==================================\n');

                data = [];

                seqPara = getSeqParameters(scenario, approach, K);
                simPara = ...
                    getSimParameters(scenario, seqPara.Nzc, SNR, numL, UE, ...
                                     scs, Tcp, min_cfo, max_cfo);

                preambles = ...
                    collectPreambles(approach, simPara.numActive, seqPara);

                [rxSignal, selectedIdx, pTable] = ...
                    preamblesFromUE(approach, preambles, simPara, seqPara);

                detected = ...
                    detectRxSignal(approach, preambles, rxSignal, simPara, seqPara);

                result = ...
                    examResult(detected, selectedIdx, pTable, seqPara.Nzc, simPara);

            end
        end
    end
end

fclose('all');