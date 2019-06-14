classdef TimeVector < matlab.mixin.CustomDisplay
    properties
        start
        delta
        stop
        
        tvec
        N
    end
    
    methods
        function tv = TimeVector(varargin)
        end
        
        function tv = alignStop(tv)
            tv.stop = tv.start + tv.delta * floor((tv.stop - tv.start) / tv.delta);
        end
        
        function tvec = get.tvec(tv)
            tvec = tv.start : tv.delta : tv.stop;
        end
        
        function N = get.N(tv)
            N = floor((tv.stop - tv.start) / tv.delta) + 1;
        end
    end
    
    methods(Static)
        function tv = from_vector(tvec)
            tv = TimeVector();
            tv.start = tvec(1);
            tv.stop = tvec(end);
            tv.delta = tvec(2) - tvec(1);
            
            tv = tv.alignStop();
        end
        
        function tv = from_start_delta_stop(start, delta, stop)
            tv = TimeVector();
            if nargin == 3
                tv.start = start;
                tv.delta = delta;
                tv.stop = stop;
            elseif nargin == 1
                tv.start = start(1);
                tv.delta = start(2);
                tv.stop = stop(3);
            end
            
            tv = tv.alignStop();
        end
        
        function tv = from_start_delta_n(start, delta, n)
             tv = TimeVector();
             tv.start = start;
             tv.delta = delta;
             tv.stop = tv.start + delta * (n-1);
        end
    end
end