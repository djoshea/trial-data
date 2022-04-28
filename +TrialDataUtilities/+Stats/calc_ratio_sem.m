function [expAoverB, semAoverB] = calc_ratio_sem(expA, semA, expB, semB, covAB)
    % see https://stats.stackexchange.com/questions/19576/variance-of-the-reciprocal-ii/19580#19580
    % var(A/B) = (E(A) / E(B))^2 (Var(A) / E(A)^2 + Var(B) / E(B)^2 - 2*Cov(A, B)/(E(A) E(B)))

    if nargin < 5
        covAB = zeros(size(A));
    end
        
    expAoverB = expA ./ expB;
        
    varA = semA.^2;
    varB = semB.^2;
    varAoverB = (expA ./ expB).^2 .* (varA ./ expA.^2 + varB ./ expB.^2 - 2*covAB ./ (expA .* expB));
    semAoverB = sqrt(varAoverB);
end