function PDP = getPDP(Nzc, inACC, zc)
    % PDP(m) calculaltion
    PDP = zeros(Nzc, 1);    
    for m = 0:Nzc-1
        PDP(m+1) = power(abs(sum(inACC .* conj(circshift(zc, -m)))/Nzc), 2);
    end
end