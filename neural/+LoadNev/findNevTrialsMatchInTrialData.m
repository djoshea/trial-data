function [qMatchInTD, tdMatchInQ] = findNevTrialsMatchInTrialData(td, Q)
    tdTrialId = td.getParam('trialId');
    tdStruct = struct('trialId', num2cell(tdTrialId));
    [qMatchInTD, tdMatchInQ] = LoadNev.matchStructArraysByField(Q, tdStruct, {'trialId'});
    
end