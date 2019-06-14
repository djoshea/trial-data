classdef ValueRange 
    properties
        low
        high
        
        closedLow logical = false;
        closedHigh logical = true;
    end
    
    methods 
        function r = AttributeValueRange(low, high)
            if nargin == 1
                if numel(low) == 2
                    r.low = low(1);
                    r.low = low(2);
                else
                    error('Invalid input to constructor');
                end
            elseif nargin == 2
                r.low = low;
                r.high = high;
            end
        end
        
        function discretize(rset)
            error('need to implement this');
        end
    end
end
        
                
                
                
        