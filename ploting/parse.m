clc;clear;close all;

rootdir = 'C:/Users/Wan/Desktop/rawData';
filelist = dir(fullfile(rootdir, '*.csv'));

% rootdir = 'C:\Users\Wan\Desktop\rawData\differentSNR_ue_01_r25_60kHz';
% filelist = dir(fullfile(rootdir, '*.csv'));

cd(rootdir);
fid = fopen('differentSNR.csv', 'w');

for i = 1:length(filelist)
    path = strcat(filelist(i).folder, '/', filelist(i).name);
    % fprintf('%s\n', path);
    T = readtable(path);
    T.Properties.VariableNames([1:5]) = {'FalseAlarm','MissDetection','avgTauI','TimingMSE','FirstAccess'};
    idx = (T.FalseAlarm > 0);
    if (~isempty(idx))
        T(idx, :) = [];
    end


    b = split(filelist(i).name, ["_"]);
    s = sprintf('%s,%s,%s,%s', b{1}, b{2},b{3}, b{4}, b{5}, b{6}, b{7});
    %, b{3}, b{4}, b{5}, b{6}, b{7} , b{8}

    % avgMSE = mean([T.TimingMSE]);
    avgACC = mean([T.FirstAccess]);
    fprintf(fid, '%s,%.4f\n', s, avgACC);
    % fprintf(fid, '%s,%.4e, %.4f\n', s, avgMSE, avgACC);
end

fclose(fid);