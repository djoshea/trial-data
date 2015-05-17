classdef ProjManual < StateSpaceProjection
    
    properties
        meanCenterBases = true;
    end
    
    methods
        function proj = ProjManual(varargin)
            proj = proj@StateSpaceProjection(varargin{:}); 
        end
        
       
        function pset = preparePsetForInference(proj, pset) 
            if proj.meanCenterBases
                pset = pset.meanSubtractBases();
            end
        end
    end

    methods
        function [decoderKbyN, encoderNbyK] = computeProjectionCoefficients(proj, pset, varargin)
            decoderKbyN = proj.decoderKbyNManual;
            encoderNbyK = proj.encoderNbyKManual;
        end
    end

end
