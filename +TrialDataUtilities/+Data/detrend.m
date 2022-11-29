function x = detrend(x, nmed)
    ntap = ceil(nmed/2);
    xf = cat(1, x(1) * ones(ntap, 1), x, x(end) * ones(ntap, 1));
    xf = smoothdata(xf, 1, 'movmedian', nmed);
    x = x - xf(ntap+1:end-ntap);
end