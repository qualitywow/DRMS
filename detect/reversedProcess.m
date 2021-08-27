function [rx, rxACC] = reversedProcess(rxSignal, cpLen, Nzc, K)
    rx = zeros(K*Nzc, 1);
    rxACC = zeros(Nzc, 1); % accumulating rx (see JSAC paper)
    rxSignal = rxSignal(cpLen+1:end);

    for v = 1:K
        temp = rxSignal((v-1)*Nzc+1:v*Nzc);
        rxFFT = fft(temp, Nzc);
        rxSeq = ifft(rxFFT, Nzc);
        rx((v-1)*Nzc+1:v*Nzc) = rxSeq;
        rxACC = rxACC + rxSeq;
    end
end