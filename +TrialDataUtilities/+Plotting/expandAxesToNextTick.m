function expandAxesToNextTick(axh, varargin)

p = inputParser;
p.addOptional('axes', 'xy', @ischar);
p.parse(varargin{:});
        
if nargin < 1
    axh = gca;
end

dir = p.Results.axes;
if strcmp(dir, 'both') || strcmp(dir, 'xy')
    axes = {'x', 'y'};
elseif strcmp(dir, 'x')
    axes = {'x'};
elseif strcmp(dir, 'y')   
    axes = {'y'};
else
    error('Parameter axes must be ''x'', ''y'', ''xy'', or ''both''');
end

for iAxis = 1:numel(axes)
    ax = axes{iAxis};
    axup = upper(ax);
    
    lim = get(axh, [axup 'Lim']);
    ticks = makecol(get(axh, [axup 'Tick']));
    
    if isempty(ticks)
        continue;
    end
    
    tickDelta = abs(ticks(2) - ticks(1));

    changed = false;
    while min(ticks) > lim(1)
        % lower tick too high, go one lower
        ticks = [min(ticks) - tickDelta; ticks]; %#ok<AGROW>
        changed = true;
    end
    
    while max(ticks) < lim(2);
        % upper tick too low, go one higher
        ticks = [ticks; max(ticks) + tickDelta]; %#ok<AGROW>
        changed = true;
    end
    
    if changed
        set(axh, [axup 'Tick'], ticks, [axup 'Lim'], [min(ticks) max(ticks)]);
    end
end

end 