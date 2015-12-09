function [optimalLambda, optimalLambdas, stats] = dpca_optimizeLambda(Xfull, ...
    trials_NbyTAbyAttrbyR, meansExcluding_NbyTAbyAttrbyR, varargin)
% Xfull is N x T x ConditionAttr
% XrandomTrials is N x T x ConditionAttr x R where R is the number of
% trials to iterate over for cross-validation
% numOfTrials is N x ConditionAttr
% 
% computes optimal regularization parameter. X is the data array. Xtrial 
% is an array storing single trials. It has one extra dimension as compared 
% with X and stores individual single trial firing rates, as opposed to the 
% trial average. numOfTrials has one dimension fewer than X and for each 
% neuron and combination of parameters (without time) specifies the number 
% of available trials in X_trial. All entries have to be larger than 1.
%
% This code assumes that time parameter is stored in the last dimension of
% X. For datasets without time, some other cross-validation needs to be
% used.
%
% [optimalLambda, optimalLambdas] = dpca_optimizeLambda(...) additionally
% returns a list of optimal lambdas found separately for each
% marginalization
%
% [...] = dpca_optimizeLambda(..., 'PARAM1',val1, 'PARAM2',val2, ...) 
% specifies optional parameter name/value pairs:
%
% 'numComps'        - how many components to use overall or in each marginalization 
%                     (default: 25)
%
% 'lambdas'         - an array of lambdas to scan
%                     (default: 1e-07 * 1.5.^[0:25])
%
% 'numRep'          - how many cross-validation iterations to perform
%                     (default: 10)
%
% 'display'         - "yes" or "no". If yes, then a figure is displayed showing
%                     reconstruction errors.
%                     (default: yes)
%
% 'combinedParams'  - cell array of cell arrays specifying 
%                     which marginalizations should be added up together,
%                     e.g. for the three-parameter case with parameters
%                           1: stimulus
%                           2: decision
%                           3: time
%                     one could use the following value:
%                     {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}}.
%
% 'filename'        - if provided, reconstruction errors and optimal lambdas will
%                     be stored in this file
%
% 'method'          - three possible ways to compute the objective to be
%                     minimized:
%                        'naive'      - reconstruction error on the test
%                                       (deprecated! do not use)
%                        'training'   - use test data to reconstruct training
%                                       data (DEFAULT)
%                        'neuronwise' - reconstruction error on the test
%                                       data computed per neuron

import(getPackageImportString);

% default input parameters
options = struct('nBasesKeep', NaN, ...
                 'nBasesPerMarginalization', NaN, ...
                 'lambdas',        1e-07 * 1.5.^(0:25), ...
                 'numRep',         10,                  ...
                 'display',        'yes',               ...
                 'combinedParams', [],                  ...
                 'filename',       [],                  ...
                 'method',         'training');

% read input parameters
optionNames = fieldnames(options);
if mod(length(varargin),2) == 1
	error('Please provide propertyName/propertyValue pairs')
end
for pair = reshape(varargin,2,[])    % pair is {propName; propValue}
	if any(strcmp(pair{1}, optionNames))
        options.(pair{1}) = pair{2};
    else
        error('%s is not a recognized parameter name', pair{1})
	end
end

assert(~isnan(options.nBasesKeep));
assert(~any(isnan(options.nBasesPerMarginalization)));

assert(size(trials_NbyTAbyAttrbyR, ndims(trials_NbyTAbyAttrbyR)) >= options.numRep, ...
    'Number of random trials provided must be greater than or equal to numRep');

% numComps = options.nBasesPerMarginalization;

% tic
% numOfTrials is N x CondAttr, need it to look like Xfull which is N x T x
% condAttr
% szNumTrials = size(numOfTrials);
% N = szNumTrials(1);
% condSize = szNumTrials(2:end);
% numOfTrials = reshape(numOfTrials, [N 1 condSize]);
% Xsum = bsxfun(@times, Xfull, numOfTrials);
%Xsum = nansum(Xtrial,5);

prog = ProgressBar(options.numRep, 'Optimizing lambda');
for rep = 1:options.numRep
    prog.update(rep);
    %fprintf(['Repetition #' num2str(rep) ' out of ' num2str(options.numRep)])
    
%     Xtest = TrialDataUtilities.DPCA.dpca_getTestTrials(Xtrial, numOfTrials);
    Xtest = TensorUtils.selectAlongDimension(trials_NbyTAbyAttrbyR, ndims(trials_NbyTAbyAttrbyR), rep);
    Xtrain = TensorUtils.selectAlongDimension(meansExcluding_NbyTAbyAttrbyR, ndims(meansExcluding_NbyTAbyAttrbyR), rep);
    
    XtestCen = bsxfun(@minus, Xtest, mean(Xtest(:,:),2));
    XtestMargs = TrialDataUtilities.DPCA.dpca_marginalize(XtestCen, ...
        'combinedParams', options.combinedParams, 'ifFlat', 'yes');
    margTestVar = nan(length(XtestMargs), 1);
    for i=1:length(XtestMargs)
        margTestVar(i) = sum(XtestMargs{i}(:).^2);
    end
    
    XtrainCen = bsxfun(@minus, Xtrain, mean(Xtrain(:,:),2));
    XtrainMargs = TrialDataUtilities.DPCA.dpca_marginalize(XtrainCen, ...
        'combinedParams', options.combinedParams, 'ifFlat', 'yes');
    margTrainVar = nan(length(XtrainMargs), 1);
    for i=1:length(XtrainMargs)
        margTrainVar(i) = sum(XtrainMargs{i}(:).^2);
    end
    
    if strcmp(options.method, 'naive') || strcmp(options.method, 'neuronwise')
        margVar_toNormalize = margTestVar;
    else
        margVar_toNormalize = margTrainVar;
    end
    
    if rep == 1
        errorsMarg = nan(numel(XtestMargs), length(options.lambdas), options.numRep);
        errors = nan(length(options.lambdas), options.numRep);
    end

    progInner = ProgressBar(length(options.lambdas), 'Testing lambda values');
    for l = 1:length(options.lambdas)
        progInner.update(l);
        
        [W,V,whichMarg] = TrialDataUtilities.DPCA.dpca(Xtrain, ...
            'nBasesKeep', options.nBasesKeep, ...
            'nBasesPerMarginalization', options.nBasesPerMarginalization, ..., ...
            'combinedParams', options.combinedParams, ...
            'lambda', options.lambdas(l));
        
        cumError = 0;
        for i=1:length(XtestMargs)
            recError = 0;
            
            if strcmp(options.method, 'naive')
                recError = sum(sum((XtestMargs{i} - V(:,whichMarg==i)*W(:,whichMarg==i)'*XtestCen(:,:)).^2));

            elseif strcmp(options.method, 'training')
                recError = sum(sum((XtrainMargs{i} - V(:,whichMarg==i)*W(:,whichMarg==i)'*XtestCen(:,:)).^2));
            
            elseif strcmp(options.method, 'neuronwise')
                % diagVW = diag(diag(V(:,whichMarg==i)*W(:,whichMarg==i)'));
                diagVW = diag(sum(V(:,whichMarg==i).*W(:,whichMarg==i), 2));
                
                recError = sum(sum((XtestMargs{i} - V(:,whichMarg==i)*W(:,whichMarg==i)'*XtestCen(:,:) ...
                    + diagVW*XtestCen(:,:)).^2));
                
                % FOR DEBUGGING: computes the same thing
                % for neur = 1:size(XtestCen,1)
                %   otherN = [1:(neur-1) (neur+1):size(XtestCen,1)];
                %   recError = recError + ...
                %       sum((XtestMargs{i}(neur,:) - V(neur,whichMarg==i)*W(otherN,whichMarg==i)'*XtestCen(otherN,:)).^2);
                % end
            end
            
            errorsMarg(i, l, rep) = recError/margVar_toNormalize(i);
            cumError = cumError + recError;
        end
        
        errors(l,rep) = cumError / sum(margVar_toNormalize);
    end
    progInner.finish();
%     fprintf('\n')
end
prog.finish();

meanError = mean(errors,2);
[~, ind] = min(meanError);
optimalLambda = options.lambdas(ind);

meanErrorMarg = mean(errorsMarg(:, :,:), 3);
[~, indm] = min(meanErrorMarg, [], 2);
optimalLambdas = options.lambdas(indm);

stats.meanErrorMarg = meanErrorMarg;
stats.rawErrorsMarg = errorsMarg;
stats.lambdas = options.lambdas;

% if ~isempty(options.filename)
%     lambdas = options.lambdas;
%     numComps = options.numComps;
%     save(options.filename, 'lambdas', 'errors', 'errorsMarg', 'optimalLambda', 'optimalLambdas', 'numComps', 'timeTaken')
% end

if strcmp(options.display, 'yes')
    figure
    
    title('Relative cross-validation errors')
    xlabel('Regularization parameter, lambda')
    ylabel('Residual variance over total test variance')
    
    hold on
    hh = patch([log(options.lambdas) fliplr(log(options.lambdas))], ...
        [min(errors,[],2)' fliplr(max(errors,[],2)')], [0 0 0]);
    set(hh, 'FaceAlpha', 0.2)
    set(hh, 'EdgeColor', 'none')
    h1 = plot(log(options.lambdas), meanError, '.-k', 'LineWidth', 2);
    plot(log(options.lambdas(ind)), meanError(ind), '.k', 'MarkerSize', 30)

    colors = lines(size(meanErrorMarg,1));
    for i=1:size(meanErrorMarg,1)
        hh = patch([log(options.lambdas) fliplr(log(options.lambdas))], ...
            [squeeze(min(errorsMarg(i,:,:),[],3)) fliplr(squeeze(max(errorsMarg(i,:,:),[],3)))], ...
            colors(i,:));
        set(hh, 'FaceAlpha', 0.2)
        set(hh, 'EdgeColor', 'none')
    end
    hh = plot(log(options.lambdas), meanErrorMarg', '.-', 'LineWidth', 1);
    for i=1:size(meanErrorMarg,1)
        plot(log(options.lambdas(indm(i))), meanErrorMarg(i,indm(i)), '.k', 'MarkerSize', 20)
    end
    
    legendText = cell(length(hh)+1, 1);
    for i = 1:length(hh)
        legendText{i} = ['Marginalization #' num2str(i)];
    end
    legendText{length(hh)+1} = 'Overall';
    legend([hh; h1], legendText, 'Location', 'East')
    
    xticks = [1e-07:1e-07:1e-06 2e-06:1e-06:1e-05 2e-05:1e-05:1e-04 2e-04:1e-04:1e-03];
    xtickLabels = num2cell(xticks);
    for i=setdiff(1:length(xticks), [1 10 19 28])
        xtickLabels{i} = '';
    end
    set(gca,'XTick', log(xticks))
    set(gca,'XTickLabel', xtickLabels)
    
    plot(xlim, [1 1], 'k')
    if numel(options.lambdas) > 1
        xlim([log(options.lambdas(1)) log(options.lambdas(end))]);
    end
    ylim([0 1.2]);
end
