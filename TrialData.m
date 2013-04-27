classdef TrialData
% TrialData represents a collection of trials, whose data is accessed via
% a TrialDataInterface. 
%
% TrialData is not a handle class, meaning that methods which modify this instance
% will return the new TrialData. Any changes made to the underlying TrialDataStore
% will be done via copy-on-write, so that changes will not propagate to another
% TrialData instance.
    
    % will eventually change to Access=protected
    properties(SetAccess=protected)
        tdi % TrialDataInterface handle: to original trial data
    end

    properties(Dependent) % Read-through to TrialDataInterface 
        meta  
        name 
        nTrials
        channels
        channelDescriptors
    end
    
    methods % get. accessors for above properties which simply refer to tdi.?
        function meta = get.meta(td)
            meta = td.tdi.meta;
        end
        
        function name = get.name(td)
            name = td.tdi.name;
        end

        function nTrials = get.nTrials(td)
            nTrials = td.tdi.nTrials;
        end

        function channels = get.channels(td)
            channels = td.tdi.channels;
        end

        function channelDescriptors = get.channelDescriptors(td)
            channelDescriptors = td.tdi.channelDescriptors;
        end
    end
    
    methods
        function ts = TrialData(tdi)
            assert(isa(tdi, 'TrialDataInterace'), 'Argument must be a TrialDataInterface instance');
                ts.tdi = tdi;
            else
        end
    end
    
end

