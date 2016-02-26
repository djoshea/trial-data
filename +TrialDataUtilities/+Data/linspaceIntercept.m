function v = linspaceIntercept(start, gap, stop, intercept)
% v = linspaceIntercept(start, gap, stop, intercept)
%
% like (start:gap:stop)', although proceed in either direction from intercept
% so that the time vector is aligned with intercept, i.e. every entry is
% equal to intercept + n*gap for integer n

v = [fliplr((intercept-gap):-gap:start), intercept:gap:stop]';
v = v(v >= start & v <= stop);


end