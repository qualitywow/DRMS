function d = getDistance(hSat, alpha)
    R = 6371e3; % m
    d = sqrt((R * sind(alpha))^2 + (hSat)^2 + 2*hSat*R) - (R * sind(alpha));
end