function solutions = SolveCongruence(a,c,m)
% Solve a linear congruence
% INPUT: a,c,m such that a*x == c(mod m)
% OUTPUT: Solutions, if any, to the congruence.
    [g,u0,v0] = EuclidMatrix(m,a);
    if ( mod(c,g) ) ~= 0
        disp('No solutions.')
        solutions = [];
        return
    end

    % disp([ num2str(m) 'x + ' num2str(a) 'y = ' num2str(g)  ])
    % disp(['Solution: (' num2str(u0) ')' num2str(m)...
    %   ' + (' num2str(v0) ')' num2str(a) ' = ' num2str(g) ]);

    u = u0 -  (a / g);
    v = v0 +  (m / g);
    x = v * (c / g);
    y = u * (c / g);
    S = @(k) x + k * (m/g);
    solutions = mod(S(0:g-1),m);
    % disp(' Solutions');
    % disp(solutions');
end

% Euclidean Algorithm, Beazout's Coefficients are stored in matrix
function [g,u,v] = EuclidMatrix(a,b)
    M = [ 1 0; 0 1 ];
    n = 0;
    while (b ~= 0)
        q = floor(a/b);
        M = M * [ q 1; 1 0];
        t = a;
        a = b;
        b = t - q * b;
        n = n + 1;
    end
    g = a;
    u = ( -1 )^n * M(2,2);
    v = (-1)^(n+1) * M(1,2);
    % disp([ num2str(u) 'x + ' num2str(v) 'y = ' num2str(g)  ])
end