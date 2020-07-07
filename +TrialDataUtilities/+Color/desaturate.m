function cmap = desaturate(cmap, amount)
    if nargin < 2
        amount = 1;
    end
    
    Kn = 18;
    lch = TrialDataUtilities.Color.convert('RGB->LCH', cmap);

    lch(:, 2) = lch(:, 2) - Kn * amount;
    lch(lch(:, 2) < 0, 2) = 0;
    
    cmap = TrialDataUtilities.Color.convert('LCH->RGB', lch);

end

% function cvec = desaturate(c, gamma)
%     if nargin < 2
%         gamma = 0.2;
%     end
%     c = TrialDataUtilities.Color.toRGB(c);
%     cvec = 1 - ((1-c) .* (1-gamma));
% end
