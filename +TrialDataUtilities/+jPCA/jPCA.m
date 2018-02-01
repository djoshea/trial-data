% This is a lightly modified version from the jPCA codepack by Mark
% Churchland and John Cunningham
%
% [Projection, Summary] = jPCA(Data, analyzeTimes, params)
%
% OUTPUTS:
%   Projection is a struct with one element per condition.
%   It contains the following fields:
%       .proj           The projection into the 6D jPCA space.
%       .times          Those times that were used
%       .projAllTimes   Pojection of all the data into the jPCs (which were derived from a subset of the data)
%       .allTimes       All the times (exactly the same as Data(1).times)
%       .tradPCAproj    The traditional PCA projection (same times as for 'proj')
%       .tradPCAprojAllTimes   Above but for all times.
%
%   Summary contains the following fields:
%       .jPCs           The jPCs in terms of the PCs (not the full-D space)
%       .PCs            The PCs (each column is a PC, each row a neuron)
%       .jPCs_highD     The jPCs, but in the original high-D space.  This is just PCs * jPCs, and is thus of size neurons x numPCs
%       .varCaptEachJPC The data variance captured by each jPC
%       .varCaptEachPC  The data variance captured by each PC
%       .R2_Mskew_2D    Fit quality (fitting dx with x) provided by Mskew, in the first 2 jPCs.
%       .R2_Mbest_2D    Same for the best M
%       .R2_Mskew_kD    Fit quality (fitting dx with x) provided by Mskew, in all the jPCs (the number of which is set by 'numPCs'.
%       .R2_Mbest_kD    Same for the best M
%
%       There are some other useful currently undocumented (but often self-explanatory) fields in
%       'Summary'
%
%   Summary also contains the following, useful for projecting new data into the jPCA space.  Also
%   useful for getting from what is in 'Projection' back into the high-D land of PSTHs
%   Note that during preprocessing we FIRST normalize, and then (when doing PCA) subtract the
%   overall mean firing rate (no relationship to the cross-condition mean subtraction).  For new
%   data you must do this in the same order (using the ORIGINAL normFactors and mean firing rats.
%   To get from low-D to high D you must do just the reverse: add back the mean FRs and the multiply
%   by the normFactors.
%       .preprocessing.normFactors = normFactors;  
%       .preprocessing.meanFReachNeuron = meanFReachNeuron; 
%
% INPUTS:
% The input 'Data' needs to be a struct, with one entry per condition.
% For a given condition, Data(c).A should hold the data (e.g. firing rates).
% Each column of A corresponds to a neuron, and each row to a timepoint.
%
% Data(c).times is an optional field.  If you provide it, only those entries that match
% 'analyzeTimes' will be used for the analysis. If  analyzeTimes == [], all times will be used.
%  
%  If you don't provide it, a '.times' field
% is created that starts at 1.  'analyzeTimes' then refers to those times.
% 
% 'params' is optional, and can contain the following fields:
%   .numPCs        Default is 6. The number of traditional PCs to use (all jPCs live within this space)
%   .normalize     Default is 'true'.  Whether or not to normalize each neurons response by its FR range.
%   .softenNorm    Default is 10.  Determines how much we undernormalize for low FR neurons.  0 means
%                  complete normalization.  10 means a neuron with FR range of 10 gets mapped to a range of 0.5.
%   .meanSubtract  Default is true.  Whether or not we remove the across-condition mean from each
%                  neurons rate.
%   .suppressBWrosettes    if present and true, the black & white rosettes are not plotted
%   .suppressHistograms    if present and true, the blue histograms are not plotted
%   .suppressText          if present and true, no text is output to the command window 
%
%  As a note on projecting more data into the same space.  This can be done with the function
%  'projectNewData'.  However, if you are going to go this route then you should turn OFF
%  'meanSubtract'.  If you wish to still mean subtract the data you should do it by hand yourself.
%  That way you can make a principled decision regarding how to treat the original data versus the
%  new-to-be-projected data.  For example, you may subtract off the mean manually across 108
%  conditions.  You might then leave out one condition and compute the jCPA plane (with meanSubtract set to false).  
%  You could then project the remaining condition into that plane using 'projectNewData'.
% 
% 
function coeffs_KbyN = jPCA(dataNbyCbyT, numPCs, varargin)

p = inputParser();
p.addParameter('suppressBWrosettes', true, @islogical);
p.addParameter('suppressHistograms', true, @islogical);
p.addParameter('suppressText', true, @islogical);
p.parse(varargin{:});

N = size(dataNbyCbyT, 1);
C = size(dataNbyCbyT, 2);
T = size(dataNbyCbyT, 3);

if rem(numPCs,2)>0
    error('you MUST ask for an even number of PCs.'); 
end

dataCTbyN = TensorUtils.reshapeByConcatenatingDims(dataNbyCbyT, {[3 2], 1});

% now do traditional PCA
[pcaMat, pcScores] = pca(dataCTbyN, 'NumComponents', numPCs, 'Economy', true, 'Centered', true);  % apply PCA to the analyzed times

% reshape back into c by t by n so we can diff
pcScoresUnpacked = reshape(pcScores, [T C numPCs]);
pcScoresFirst = reshape(pcScoresUnpacked(1, :, :), [C numPCs]);
diffTensor = diff(pcScoresUnpacked, 1, 1);
preTensor = pcScoresUnpacked(1:end-1, :, :);

% and back into CT by N
dState = reshape(diffTensor, [C*(T-1) numPCs]);
preState = reshape(preTensor, [C*(T-1) numPCs]);
tMaskValid = all(~isnan(dState) & ~isnan(preState), 2);
dState = dState(tMaskValid, :);
preState = preState(tMaskValid, :);

% GET M & Mskew
% compute dState, and use that to find the best M and Mskew that predict dState from the state
% these are used to take the derivative

% we are interested in the eqn that explains the derivative as a function of the state: dState/dt = M*State

% first compute the best M (of any type)
% note, we have converted dState and Ared to have time running horizontally
M = (dState'/preState');  % M takes the state and provides a fit to dState
% Note on sizes of matrices:
% dState' and preState' have time running horizontally and state dimension running vertically 
% We are thus solving for dx = Mx.
% M is a matrix that takes a column state vector and gives the derivative

% now compute Mskew using John's method
% Mskew expects time to run vertically, transpose result so Mskew in the same format as M
% (that is, Mskew will transform a column state vector into dx)
Mskew = TrialDataUtilities.jPCA.skewSymRegress(dState,preState)';  % this is the best Mskew for the same equation

%% USE Mskew to get the jPCs

% get the eigenvalues and eigenvectors
[V,D] = eig(Mskew); % V are the eigenvectors, D contains the eigenvalues
evals = diag(D); % eigenvalues

% the eigenvalues are usually in order, but not always.  We want the biggest
[~,sortIndices] = sort(abs(evals),1,'descend');
evals = evals(sortIndices);  % reorder the eigenvalues
evals = imag(evals);  % get rid of any tiny real part
V = V(:,sortIndices);  % reorder the eigenvectors (base on eigenvalue size)

% Eigenvalues will be displayed to confirm that everything is working
% unless we are asked not to output text
if ~exist('params', 'var') || ~isfield(params,'suppressText') || ~params.suppressText
    disp('eigenvalues of Mskew: ');
    for i = 1:length(evals)
        if evals(i) > 0
            fprintf('                  %1.3fi', evals(i));
        else
            fprintf('     %1.3fi \n', evals(i));
        end
    end
end

jPCs = zeros(size(V));
for pair = 1:numPCs/2
    vi1 = 1+2*(pair-1);
    vi2 = 2*pair;
    
    VconjPair = V(:,[vi1,vi2]);  % a conjugate pair of eigenvectors
    evConjPair = evals([vi1,vi2]); % and their eigenvalues
    VconjPair = getRealVs(VconjPair,evConjPair, pcScoresFirst, pcScores);
    
    jPCs(:,[vi1,vi2]) = VconjPair;
end

% pcaMat is K by N
% jPCS is K by K
% data is 
coeffs_KbyN = (pcaMat * jPCs)'; 

end
% 
% %% Get the projections
% 

% TxK * K*K == T*N * (N*K * K*K)
% output is K*N

% proj = Ared * jPCs;
% projAllTimes = bigAred * jPCs;
% tradPCA_AllTimes = bsxfun(@minus, bigA, mean(smallA)) * PCvectors;  % mean center in exactly the same way as for the shorter time period.
% crossCondMeanAllTimes = meanAred * jPCs;
% 
% % Do some annoying output formatting.
% % Put things back so we have one entry per condition
% index1 = 1;
% index2 = 1;
% for c = 1:numConds
%     index1b = index1 + numAnalyzedTimes -1;  % we will go from index1 to this point
%     index2b = index2 + numTimes -1;  % we will go from index2 to this point
%     
%     Projection(c).proj = proj(index1:index1b,:);
%     Projection(c).times = Data(1).times(analyzeIndices);
%     Projection(c).projAllTimes = projAllTimes(index2:index2b,:);
%     Projection(c).allTimes = Data(1).times;
%     Projection(c).tradPCAproj = Ared(index1:index1b,:);
%     Projection(c).tradPCAprojAllTimes = tradPCA_AllTimes(index2:index2b,:);
%     
%     index1 = index1+numAnalyzedTimes;
%     index2 = index2+numTimes;
% end
%    
% %% Done computing the projections, plot the rosette
% 
% % do this unless params contains a field 'suppressBWrosettes' that is true
% if ~exist('params', 'var') || ~isfield(params,'suppressBWrosettes') || ~params.suppressBWrosettes
%     plotRosette(Projection, 1);  % primary plane
%     plotRosette(Projection, 2);  % secondary plane
% end
% 
% 
% %% SUMMARY STATS
% %% compute R2 for the fit provided by M and Mskew
% 
% % R2 Full-D
% fitErrorM = dState'- M*preState';
% fitErrorMskew = dState'- Mskew*preState';
% varDState = sum(dState(:).^2);  % original data variance
% 
% R2_Mbest_kD = (varDState - sum(fitErrorM(:).^2)) / varDState;  % how much is explained by the overall fit via M
% R2_Mskew_kD = (varDState - sum(fitErrorMskew(:).^2)) / varDState;  % how much by is explained via Mskew
% 
% % unless asked to not output text
% if ~exist('params', 'var') || ~isfield(params,'suppressText') || ~params.suppressText
%     fprintf('%% R^2 for Mbest (all %d dims):   %1.2f\n', numPCs, R2_Mbest_kD);
%     fprintf('%% R^2 for Mskew (all %d dims):   %1.2f  <<---------------\n', numPCs, R2_Mskew_kD);
% end
% 
% 
% % R2 2-D primary jPCA plane
% fitErrorM_2D = jPCs(:,1:2)' * fitErrorM;  % error projected into the primary plane
% fitErrorMskew_2D = jPCs(:,1:2)' * fitErrorMskew;  % error projected into the primary plane
% dState_2D = jPCs(:,1:2)' * dState'; % project dState into the primary plane
% varDState_2D = sum(dState_2D(:).^2); % and get its variance
% 
% R2_Mbest_2D = (varDState_2D - sum(fitErrorM_2D(:).^2)) / varDState_2D;  % how much is explained by the overall fit via M
% R2_Mskew_2D = (varDState_2D - sum(fitErrorMskew_2D(:).^2)) / varDState_2D;  % how much by is explained via Mskew
% 
% if ~exist('params', 'var') || ~isfield(params,'suppressText') || ~params.suppressText
%     fprintf('%% R^2 for Mbest (primary 2D plane):   %1.2f\n', R2_Mbest_2D);
%     fprintf('%% R^2 for Mskew (primary 2D plane):   %1.2f  <<---------------\n', R2_Mskew_2D);
% end
% 
% %% variance catpured by the jPCs
% origVar = sum(sum( bsxfun(@minus, smallA, mean(smallA)).^2));
% varCaptEachPC = sum(Ared.^2) / origVar;  % this equals latent(1:numPCs) / sum(latent)
% varCaptEachJPC = sum((Ared*jPCs).^2) / origVar;
% varCaptEachPlane = reshape(varCaptEachJPC, 2, numPCs/2);
% varCaptEachPlane = sum(varCaptEachPlane);
% 
% 
% %% Analysis of whether things really look like rotations (makes plots)
% 
% for jPCplane = 1:2
%     phaseData = getPhase(Projection, jPCplane);  % does the key analysis
%     
%     if exist('params', 'var')
%         cstats = plotPhaseDiff(phaseData, params, jPCplane);  % plots the histogram.  'params' is just what the user passed, so plots can be suppressed
%     else
%         cstats = plotPhaseDiff(phaseData, [], jPCplane);
%     end
%     
%     if jPCplane == 1
%         circStats = cstats;  % keep only for the primary plane
%     end
% end
% 
% %% Make the summary output structure
% 
% Summary.jPCs = jPCs;
% Summary.PCs = PCvectors;
% Summary.jPCs_highD = PCvectors * jPCs;
% Summary.varCaptEachJPC = varCaptEachJPC;
% Summary.varCaptEachPC = varCaptEachPC;
% Summary.varCaptEachPlane = varCaptEachPlane;
% Summary.Mbest = M;
% Summary.Mskew = Mskew;
% Summary.fitErrorM = fitErrorM;
% Summary.fitErrorMskew = fitErrorMskew;
% Summary.R2_Mskew_2D = R2_Mskew_2D;
% Summary.R2_Mbest_2D = R2_Mbest_2D;
% Summary.R2_Mskew_kD = R2_Mskew_kD;
% Summary.R2_Mbest_kD = R2_Mbest_kD;
% Summary.circStats = circStats;
% Summary.acrossCondMeanRemoved = meanSubtract;
% Summary.crossCondMean = crossCondMeanAllTimes(analyzeIndices,:);
% Summary.crossCondMeanAllTimes = crossCondMeanAllTimes;
% Summary.preprocessing.normFactors = normFactors;  % Used for projecting new data from the same neurons into the jPC space
% Summary.preprocessing.meanFReachNeuron = meanFReachNeuron; % You should first normalize and then mean subtract using this (the original) mean
% % conversely, to come back out, you must add the mean back on and then MULTIPLY by the normFactors
% % to undo the normalization.
% 
% %% DONE (only functions below)
% 
% 
% %% Inline function that gets the real analogue of the eigenvectors
function Vr = getRealVs(V, evals, pcScoresFirst, pcScores)

    % get real vectors made from the eigenvectors
    
    % by paying attention to this order, things will always rotate CCW
    if abs(evals(1))>0  % if the eigenvalue with negative imaginary component comes first
        Vr = [V(:,1) + V(:,2), (V(:,1) - V(:,2))*1i]; 
    else
        Vr = [V(:,2) + V(:,1), (V(:,2) - V(:,1))*1i];
    end
    Vr = Vr / sqrt(2);

    % now get axes aligned so that plan is spread mostly along the horizontal axis
    testProj = (Vr'*pcScoresFirst')'; % just picks out the plan times
    rotV = pca(testProj);
    crossProd = cross([rotV(:,1);0], [rotV(:,2);0]);
    if crossProd(3) < 0, rotV(:,2) = -rotV(:,2); end   % make sure the second vector is 90 degrees clockwise from the first
    Vr = Vr*rotV; 

    % flip both axes if necessary so that the maximum move excursion is in the positive direction
    testProj = (Vr'*pcScores')';  % all the times
    if max(abs(testProj(:,2))) > max(testProj(:,2))  % 2nd column is the putative 'muscle potent' direction.
        Vr = -Vr;
    end
end
% 
% end
% %% end of main function
% 
% 
% %% Some inline functions are below
% 
% %% Plot the rosette itself (may need to spruce this up and move to a subfunction)
% function plotRosette(Proj, whichPair)
% 
%     d1 = 1 + 2*(whichPair-1);
%     d2 = d1+1;
% 
%     numConds = length(Proj);
% 
%     figure;
% 
%     % first deal with the ellipse for the plan variance (we want this under the rest of the data)
%     planData = zeros(numConds,2);
%     for c = 1:numConds
%         planData(c,:) = Proj(c).proj(1,[d1,d2]);
%     end
%     planVars = var(planData);
%     circle([0 0], 2*planVars.^0.5, 0.6*[1 1 1], 1); hold on;
%     %fprintf('ratio of plan variances = %1.3f (hor var / vert var)\n', planVars(1)/planVars(2));
% 
%     allD = vertcat(Proj(:).proj);  % just for getting axes
%     allD = allD(:,d1:d2);
%     mxVal = max(abs(allD(:)));
%     axLim = mxVal*1.05*[-1 1 -1 1];
%     arrowSize = 5;
%     for c = 1:numConds
%         plot(Proj(c).proj(:,d1), Proj(c).proj(:,d2), 'k');
%         plot(Proj(c).proj(1,d1), Proj(c).proj(1,d2), 'ko', 'markerFaceColor', [0.7 0.9 0.9]);
% 
%         penultimatePoint = [Proj(c).proj(end-1,d1), Proj(c).proj(end-1,d2)];
%         lastPoint = [Proj(c).proj(end,d1), Proj(c).proj(end,d2)];
%         arrowMMC(penultimatePoint, lastPoint, [], arrowSize, axLim);
% 
%     end
% 
%     axis(axLim);
%     axis square;
%     plot(0,0,'k+');
%     
%     title(sprintf('jPCA plane %d', whichPair));
% end
% 
% 
% %% Getting the phases
% function phaseData = getPhase(Proj, whichPair)
%     numConds = length(Proj);
%     d1 = 1 + 2*(whichPair-1);
%     d2 = d1+1;
%     
%     for c=1:numConds
%         data = Proj(c).proj(:,[d1,d2]);
%         phase = atan2(data(:,2), data(:,1));  % Y comes first for atan2
%         
%         deltaData = diff(data);
%         phaseOfDelta = atan2(deltaData(:,2), deltaData(:,1));  % Y comes first for atan2
%         phaseOfDelta = [phaseOfDelta(1); phaseOfDelta];  %#ok<AGROW> % so same length as phase
%         radius = sum(data.^2,2).^0.5;
%         
%         % collect and format
%         % make things run horizontally so they can be easily concatenated.
%         phaseData(c).phase = phase'; %#ok<AGROW>
%         phaseData(c).phaseOfDelta = phaseOfDelta'; %#ok<AGROW>
%         phaseData(c).radius = radius'; %#ok<AGROW>
%         
%         % angle between state vector and Dstate vector
%         % between -pi and pi
%         phaseData(c).phaseDiff = minusPi2Pi(phaseData(c).phaseOfDelta - phaseData(c).phase); %#ok<AGROW>
%     end
%     
% end
% 
% 
% 
% %% plotting the phase difference between dx(t)/dt and x(t) where x is the 2D state
% function circStatsSummary = plotPhaseDiff(phaseData, params, jPCplane)
%     % compute the circular mean of the data, weighted by the r's
%     circMn = circ_mean([phaseData.phaseDiff]', [phaseData.radius]');
%     resultantVect = circ_r([phaseData.phaseDiff]', [phaseData.radius]');
%     
%     
%     bins = pi*(-1:0.1:1);
%     cnts = histc([phaseData.phaseDiff], bins);  % not for plotting, but for passing back out
%     
% 
%     % do this unless params contains a field 'suppressHistograms' that is true
%     if ~exist('params', 'var') || ~isfield(params,'suppressHistograms') || ~params.suppressHistograms
%         figure;
%         hist([phaseData.phaseDiff], bins); hold on;
%         plot(circMn, 20, 'ro', 'markerFa', 'r', 'markerSiz', 8);
%         plot(pi/2*[-1 1], [0 0], 'ko', 'markerFa', 'r', 'markerSiz', 8);
%         set(gca,'XLim',pi*[-1 1]);
%         title(sprintf('jPCs plane %d', jPCplane));
%     end
%     
%     %fprintf('(pi/2 is %1.2f) The circular mean (weighted) is %1.2f\n', pi/2, circMn);
%     
%     % compute the average dot product of each datum (the angle difference for one time and condition)
%     % with pi/2.  Will be one for perfect rotations, and zero for random data or expansions /
%     % contractions.
%     avgDP = averageDotProduct([phaseData.phaseDiff]', pi/2);
%     %fprintf('the average dot product with pi/2 is %1.4f  <<---------------\n', avgDP);
%     
%     circStatsSummary.circMn = circMn;
%     circStatsSummary.resultantVect = resultantVect;
%     circStatsSummary.avgDPwithPiOver2 = avgDP;  % note this basically cant be <0 and definitely cant be >1
%     circStatsSummary.DISTRIBUTION.bins = bins;
%     circStatsSummary.DISTRIBUTION.binCenters = (bins(1:end-1) + bins(2:end))/2;
%     circStatsSummary.DISTRIBUTION.counts = cnts(1:end-1);
%     circStatsSummary.RAW.rawData = [phaseData.phaseDiff]';
%     circStatsSummary.RAW.rawRadii = [phaseData.radius]';
%     
% end
% 
% %% for computing the average dot product with a comparison angle
% function avgDP = averageDotProduct(angles, compAngle, varargin)
%     
%     x = cos(angles-compAngle);
%     
%     if ~isempty(varargin)
%         avgDP = mean(x.*varargin{1}) / mean(varargin{1});  % weighted sum
%     else
%         avgDP = mean(x);
%     end
% end
%         
% 















