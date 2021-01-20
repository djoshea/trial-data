classdef Run < LowRankLFADS.Run
    properties
        trialDataSet % loaded from DatasetCollection
        gpfaSequenceData % loaded from GPFA run
    end
    
    properties(Dependent)
        pathGPFAOutput
    end
    
    % You must provide implementation for these
    methods(Abstract) 
        td = prepareTrialDataForLFADS(r, td, varargin) 
        out = generateCountsForTrialData(r, ds, trialData, varargin);
        out = generateCountsForTrialDataForAlignment(r, ds, trialData, varargin); % this must return trials of duration
        names = getNeuralChannelNamesFromTrialData(r, ds, td);
    end
    
    % Override as necessary
    methods 
        function td = prepareTrialDataForAlignment(r, td)
            % if you override this, set usesDifferentTrialDataForAlignment
            % to return true or this method wont be called
            td = r.prepareTrialDataForLFADS(td);
        end
        
        function chList = listChannelsForLFADS(r, td, varargin) %#ok<INUSD>
            warning('Override listChannelsForLFADS in order to specify channel names in the posterior mean rates');
            chList = {};
        end
    end      
    
    methods % needed for LFADS prep
        function r = Run(varargin) 
           r@LowRankLFADS.Run(varargin{:});
        end

        % counts : nTrials x nChannels x nTime tensor
        function trials = generateCountsForDataset(r, dataset, varargin)
            td = dataset.loadData();
            
            td = r.prepareTrialDataForLFADS(td);
            
            if ~td.hasParamChannel('trialId')
                warning('Provide parameter channel trialId to uniquely identify each trial');
                td = td.addParamChannel('trialId', (1:td.nTrials)');
            end
            td = td.selectValidTrials();
            
            trials = r.generateCountsForTrialData(dataset, td, varargin);
            trials = makecol(trials);
        end
        
        function trials = generateCountsForDatasetForAlignment(r, dataset, varargin)
            td = dataset.loadData();
            
            td = r.prepareTrialDataForAlignment(td);
            
            if ~td.hasParamChannel('trialId')
                warning('Provide parameter channel trialId to uniquely identify each trial');
                td = td.addParamChannel('trialId', (1:td.nTrials)');
            end
            td = td.selectValidTrials();
            
            trials = r.generateCountsForTrialDataForAlignment(dataset, td, varargin);
            trials = makecol(trials);
        end
        
        function names = getNeuralChannelNamesForDataset(r, dataset)
            td = dataset.loadData();
            names = r.getNeuralChannelNamesFromTrialData(dataset, td);
        end   
    end
    
    methods % Working with LFADS generated data with TrialData
        function tdSet = loadTrialDataFromDatasetCollection(r, varargin)
            p = inputParser();
            p.addOptional('reload', false, @islogical);
            p.addParameter('datasetIdx', 1:r.nDatasets, @isvector);
            p.parse(varargin{:});
            datasetIdx = LowRankLFADS.Utils.vectorMaskToIndices(p.Results.datasetIdx);
            
            if ~isempty(r.trialDataSet) && ~p.Results.reload
                tdSet = r.trialDataSet(datasetIdx);
                return;
            end
            
            tdSet = cell(numel(datasetIdx), 1);
            
            prog = ProgressBar(numel(datasetIdx), 'Loading trialData from datasets');
            for iiDS = 1:numel(datasetIdx)
                tdSet{iiDS} = r.datasets(datasetIdx(iiDS)).loadData();
                
            end
            prog.finish();
            
            % set this here as prepareTrialDataForLFADS may need to see the full set
            r.trialDataSet = tdSet;
            
            prog = ProgressBar(numel(datasetIdx), 'Preparing trialData');
            for iiDS = 1:numel(datasetIdx)
                tdSet{iiDS} = r.prepareTrialDataForLFADS(tdSet{iiDS});
                prog.update(iiDS);
            end
            prog.finish();
            
            if isequal(datasetIdx, (1:r.nDatasets)')
                r.trialDataSet = tdSet;
            end
        end
        
        function tdSet = addPosteriorMeansToTrialData(r, varargin)
            p = inputParser();
            p.addParameter('datasetIdx', 1:r.nDatasets, @isvector);
            p.addParameter('addRates', true, @islogical);
            p.addParameter('addControllerOutputs', true, @islogical);
            p.addParameter('addGeneratorICs', true, @islogical); 
            p.parse(varargin{:});
            datasetIdx = LowRankLFADS.Utils.vectorMaskToIndices(p.Results.datasetIdx);
            
            trialDataSet = r.loadTrialDataFromDatasetCollection('datasetIdx', datasetIdx); %#ok<*PROPLC> % these will be prepared for LFADS
            [pms, pms_valid] = r.loadPosteriorMeans('datasetIdx', datasetIdx);
            
            tdSet = cellvec(numel(datasetIdx));
            if ~all(pms_valid)
                warning('Posterior means not found for %d runs', nnz(~pms_valid));
                return;
            end
            
            timeField = 'posteriorMeans_time';
            prog = ProgressBar(numel(datasetIdx), 'Merging posterior mean data into trialData for each dataset');
            for iiDS = 1:numel(datasetIdx)
                prog.update(iiDS);
                     
                td = trialDataSet{iiDS};
                inflate = @(data) TensorUtils.inflateMaskedTensor(data, 1, td.valid); % loaded data only spans valid trials in td.data
                inflate2 = @(data) TensorUtils.inflateMaskedTensor(data, 2, td.valid);
                pm = pms(iiDS);
                
                % data must be nTrials x nTime x nChannels tensor
                td = td.dropChannels({'controllerOutputs', 'factors', 'generatorStates', 'rates'});
                td = td.dropChannels({'generatorIC', 'post_g0_mean', 'post_g0_logvar'});
                
                facNames = {}; % genNames('f', pm.nFactors)
                td = td.addAnalogChannelGroup('factors', inflate(permute(pm.factors, [3 2 1])), pm.time, ...
                    'subChannelNames', facNames', 'timeField', timeField, 'isAligned', true);
                    
                genStateNames = {}; % genNames('g', pm.nGeneratorUnits)
                td = td.addAnalogChannelGroup('generatorStates', inflate(permute(pm.generator_states, [3 2 1])), pm.time, ...
                    'subChannelNames', genStateNames, 'timeField', timeField, 'isAligned', true);
                
                rnames = r.listChannelsForLFADS(td);
                if isempty(rnames) || numel(rnames) == 1
                    rnames = {};
                else
                    rnames = cellfun(@(n) sprintf('rate_%s', n), rnames, 'UniformOutput', false);
                end
                
                if p.Results.addRates
                    td = td.addAnalogChannelGroup('rates', inflate(permute(pm.rates, [3 2 1])), pm.time, 'timeField', timeField, ...
                        'subChannelNames', rnames, 'units', 'spikes/sec', 'isAligned', true);
                end
                
                if pm.nControllerOutputs > 0 && p.Results.addControllerOutputs
                    coNames = {}; % genNames('co', pm.nControllerOutputs)
                    td = td.addAnalogChannelGroup('controllerOutputs', inflate(permute(pm.controller_outputs, [3 2 1])), pm.time, ...
                        'subChannelNames', coNames, 'timeField', timeField, 'isAligned', true);
                end
                
                if p.Results.addGeneratorICs
                    td = td.addVectorParamAccessAsMatrix('generatorIC', TensorUtils.splitAlongDimension(inflate2(pm.generator_ics), 2)');
                end
                td = td.addVectorParamAccessAsMatrix('post_g0_mean', TensorUtils.splitAlongDimension(inflate2(pm.post_g0_mean), 2)');
                td = td.addVectorParamAccessAsMatrix('post_g0_logvar', TensorUtils.splitAlongDimension(inflate2(pm.post_g0_logvar), 2)');
                
                tdSet{iiDS} = td;
            end
            
            if isequal(datasetIdx, (1:r.nDatasets)')
                r.trialDataSet = tdSet;
            end
        end
        
        function exportTrainedTrialData(r, varargin)
            p = inputParser();
            p.addParameter('exportPath', fullfile(r.path, 'export_trialData'), @ischar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            exportPath = p.Results.exportPath;
            if ~exist(exportPath, 'dir')
                mkdir(exportPath);
            end
            
            fprintf('Exporting to %s\n', exportPath);
            
            mtp = r.loadModelTrainedParams();
            mtp.exportToHDF5(fullfile(exportPath, 'modelTrainedParams.h5'));
            
            prog = LowRankLFADS.Utils.ProgressBar(r.nDatasets, 'Exporting trial data with posterior means');
            for iDS = 1:r.nDatasets
                prog.update(iDS);
                dsname = r.datasetNames{iDS};
                td_path = fullfile(exportPath, sprintf('td_%s', dsname));
                
                tds = r.addPosteriorMeansToTrialData('datasetIdx', iDS, p.Unmatched);
                tds{1}.saveFast(td_path);
            end
            prog.finish();
        end

        function tdSet = getTrialDataWithVirtualPopulationRate(r, dsIdx, W, b)
            timeField = 'posteriorMeans_time';
            
            % concatenated readout matrices for each dataset
            ro = r.loadReadoutMatricesByDataset();
            
            if ~exist('W', 'var')
                W = cat(1, ro.rates_W); % N_all x Far')
                b = cat(1, ro.rates_b); % N_all x 1
            end
            
            % row normalize as is done in lfads code
            Wnorm = W ./ sqrt(sum(W.^2, 2));
            N_all = size(Wnorm, 1);
            
            vrNames = genNames('vr', N_all);
            dsIdx = TensorUtils.vectorMaskToIndices(dsIdx);
            
            prog = ProgressBar(numel(dsIdx), 'Adding virtual population rates to trialData for each dataset');
            tdSet = cellvec(numel(dsIdx));
            for i = 1:numel(tdSet)
                idx = dsIdx(i);
                prog.update(idx);
                td =  r.trialDataSet{idx};
                pm = r.posteriorMeans(idx);
                
                % data must be nTrials x nTime x nChannels tensor
                td = td.dropAnalogChannelGroup({'virtualRates'});
                
                virtualRates = exp(bsxfun(@plus, TensorUtils.linearCombinationAlongDimension(pm.factors, 1, Wnorm), b));
                
                td = td.addAnalogChannelGroup('virtualRates', vrNames, ...
                    permute(virtualRates, [3 2 1]), pm.time, 'timeField', timeField, 'isAligned', true);
                tdSet{i} = td;
            end
            
            % automatically align trial data the same way as the LFADS data was
            % prepared
            tdSet = r.setupTrialDataAsLFADS(tdSet);
            
            function names = genNames(pre, n)
                names = arrayfun(@(i) sprintf('%s%i', pre, i), (1:n)', 'UniformOutput', false);
            end
        end
            
        function Tset = buildTStructs(r, varargin)
            % opts can specify :
            %   opts.lag - how much to lag the neural and kinematic data
            %     (each trial will be trimmed by lag milliseconds)
            % tStart (scalar) : where to start relative to r.params.align,
            % if NaN, use the first output from LFADS
            %   opts.neuralBinSizeMS - what is the neural data currently binned at?
            %   opts.neuralFieldName - specify the field name of the neural
            %     data. Default is 'y'

            p = inputParser();
            p.addParameter('datasetMask',truevec(r.nDatasets), @(x) true);
            p.addParameter('source', 'neural', @(x) ismember(x, {'neural', 'smoothed_neural', 'rates', 'virtualRates', ...
                'factors', 'generatorStates', 'gpfa_xorth'}));
            p.addParameter('behavioralChannels', {'handX', 'handY', 'handVelocityX', 'handVelocityY'}, @iscellstr);
            p.addParameter('neural_smooth', 40, @isscalar); % for smoothed_neural, SD of the Gaussian
            p.addParameter('rates_W', [], @ismatrix); % for virtualRates, defaults to full set of virtual neurons
            p.addParameter('rates_b', [], @isvector); % for virtualRates, defaults to full set of virtual neurons
            p.addParameter('align', 'GoCue', @ischar);
            p.addParameter('binWidth', 20, @isscalar);
            p.addParameter('tStart', NaN, @isscalar);
            p.addParameter('tStop', NaN, @isscalar);
            p.addParameter('lag', 0, @isscalar);
            p.addParameter('neuralFieldName', 'y', @ischar);
            p.addParameter('neuralBinSize', 1, @isscalar);
            p.parse(varargin{:});
            
            tLag = p.Results.lag;
            binWidth = p.Results.binWidth;
            source = p.Results.source;
            align = p.Results.align;
            
            dsIdx = TensorUtils.vectorMaskToIndices(p.Results.datasetMask);
            
            % take the limits of the posterior means in the appropriate
            % alignment
            prog = ProgressBar(r.nDatasets, 'Building T structs');
            Tset = cell(numel(dsIdx), 1);
            for iDS = 1:numel(dsIdx)
                prog.update(iDS);
                td = r.trialDataSet{dsIdx(iDS)};
            
                [~, tvec] = td.getAnalogChannelGroupAsTensor('rates', 'minTrialFraction', 1);
                tStart = nanmax(p.Results.tStart, tvec(1));
                tStop = nanmin(p.Results.tStop, tvec(end));
                
                % round to nearest binWidth multiple
                tStart = ceil(tStart / binWidth) * binWidth;
                tStop = floor(tStop / binWidth) * binWidth;
                
                % align trial data
                td = td.unalign.start(align, tStart).stop(align, tStop);
                
                % lag kinematics relative to neural data 
                % Trials x Time x 4 channels --> 4 x Time x Trials
                [kinData, time] = td.lag(tLag).getAnalogMultiAsTensor(...
                    p.Results.behavioralChannels, 'timeDelta', binWidth);
                if any(isnan(kinData(:)))
                    warning('%d nans in behavioral data', nnz(isnan(kinData(:))));
                end
                kinematics = permute(kinData, [3 2 1]);
                
                % append ones as lasst channel
                kinematics = TensorUtils.expandAlongDims(kinematics, 1, 1, 1);
                
                % split by trials
                kinematics = squeeze(TensorUtils.splitAlongDimension(kinematics, 3));
                
                switch source
                    case 'neural'
                        % Trials x Time x Neurons --> Neurons x Time x Trials
                        decode = permute(td.getSpikeBinnedCounts(td.listSpikeChannels, 'binWidthMs', binWidth), [3 2 1]);
                    
                    case 'smoothed_neural'
                        % Trials x Time x Neurons --> Neurons x Time x Trials
                        sf = GaussianSpikeFilter('sigma', p.Results.neural_smooth);
                        sf.timeDelta = binWidth;
                        decode = permute(td.getSpikeRateFilteredAsMatrix(td.listSpikeChannels, 'spikeFilter', sf), [3 2 1]);
                        mask = isnan(decode);
                        if any(mask(:))
                            warning('Setting %d values from NaN to 0', nnz(mask));
                            decode(mask) = 0;
                        end
                            
                    case 'rates'
                        % Trials x Time x Neurons --> Neurons x Time x Trials
                        decode = permute(td.getAnalogChannelGroupAsTensor('rates', 'timeDelta', binWidth), [3 2 1]);
                    
                    case 'virtualRates'
                        % Trials x Time x NeuronsAll --> Neurons x Time x Trials
                        factorTensor = permute(td.getAnalogChannelGroupAsTensor('factors', 'timeDelta', binWidth), [3 2 1]);
                        
                        % concatenated readout matrices for each dataset
                        ro = r.loadReadoutMatricesByDataset();
                        if isempty(p.Results.rates_W)
                            W = cat(1, ro.rates_W); % N_all x Factors
                            b = cat(1, ro.rates_b); % N_all x 1
                        else
                            W = p.Results.rates_W;
                            b = p.Results.rates_b;
                        end
                        Wnorm = W ./ sqrt(sum(W.^2, 2)); % row normalize as is done in lfads code  
            
                        % rates = exp(W*factors + b)
                        decode = exp(TensorUtils.linearCombinationAlongDimension(factorTensor, 1, Wnorm) + b);
                    
                    case 'factors'
                        % Trials x Time x Neurons --> Neurons x Time x Trials
                        decode = permute(td.getAnalogChannelGroupAsTensor('factors', 'timeDelta', binWidth), [3 2 1]);
                    
                    case 'generatorStates'
                        % Trials x Time x Neurons --> Neurons x Time x Trials
                        decode = permute(td.getAnalogChannelGroupAsTensor('generatorStates', 'timeDelta', binWidth), [3 2 1]);
                    
                    case 'gpfa_xorth'
                        % Trials x Time x nGPFA --> nGPFA x Time x Trials
                        decode = permute(td.getAnalogChannelGroupAsTensor('gpfa_xorth', 'timeDelta', binWidth), [3 2 1]);
                end
                
                % split by trials
                decode = squeeze(TensorUtils.splitAlongDimension(decode, 3));
                
                % number of timepoints
                nTimepoints = cellfun(@(x) size(x, 2), decode, 'UniformOutput', false);
                
                valid = num2cell(td.valid);
                
                conditions = num2cell(td.conditionIdx);
                
                % split by trials and assign into struct
                Tset{iDS} = struct('valid', valid, 'X', kinematics, 'Z', decode, 'time', time, 'T', nTimepoints, 'dt', p.Results.binWidth, 'condition', conditions);
            end
            
            prog.finish();
        end
        
        function tdSet = addDecodeTStructToTrialData(r, Tset, varargin)
            p = inputParser();
            p.addParameter('field', 'xk', @ischar);
            p.addParameter('channelPrefix', '', @ischar);
            p.addParameter('channelNames', {'positionX', 'positionY', 'velocityX', 'velocityY'}, @ischar);
            p.parse(varargin{:});
            
            % uses field xk to produce x and y velocity channels
            if isempty(r.trialDataSet)
                r.loadTrialDataFromDatasetCollection();
            end
            
            if ~iscell(Tset)
                assert(r.nDatasets == 1);
                Tset = {Tset};
            end
            
            channelPrefix = p.Results.channelPrefix;
            timeField = sprintf('%s_time', channelPrefix);
            channelNames = cellfun(@(x) sprintf('%s%s', channelPrefix, x), p.Results.channelNames, 'UniformOutput', false);
            
            prog = ProgressBar(r.nDatasets, 'Merging Tstructs into trialData for each dataset');
            tdSet = cellvec(r.nDatasets);
            for i = 1:r.nDatasets
                prog.update(i);
                td =  r.trialDataSet{i};
                T = Tset{i};
                
                % data must be nTrials x 1 cell of nTime x nChannels
                % matrices
                dataThis = arrayfun(@(s) s.xk(1:4, :)', T, 'UniformOutput', false);
                timeThis = arrayfun(@(s) s.time, T, 'UniformOutput', false);
                td = td.dropAnalogChannelGroup(channelPrefix);
                
                td = td.addAnalogChannelGroup(channelPrefix, channelNames, ...
                    dataThis, timeThis, 'timeField', timeField, 'isAligned', true);

                tdSet{i} = td;
            end
            
            % automatically align trial data the same way as the LFADS data was
            % prepared
            r.trialDataSet = tdSet;
        end
    end
    
    methods % GPFA smoothing of neural data
        function p = get.pathGPFAOutput(r)
            if isempty(r.runCollection)
                p = '';
            else
                p = fullfile(r.path, 'gpfaOutput');
            end
        end

        function p = getGpfaResultsFile(r, dtMS, nLatents)
            gpfaResultsDir = r.getGpfaResultsDir(dtMS, nLatents);
            p = fullfile(gpfaResultsDir, 'all_results.mat');
        end

        function gpfaResultsDir = getGpfaResultsDir(r, dtMS, nLatents)
            gpfaResultsDir = fullfile(r.pathGPFAOutput, ...
                                      sprintf('nLat_%03i_binSizeMS_%03i', nLatents, dtMS));
        end
        
        function resultsOut = loadGPFA(r, dtMS, nLatents)
            gpfaOutputFile = r.getGpfaResultsFile(dtMS, nLatents);
            tmp = load(gpfaOutputFile);
            resultsOut = tmp.resultsOut;
            r.gpfaSequenceData = resultsOut;
        end
        
        function resultsOut = runGPFA(r, dtMS, nLatents, deleteExistingResults)
        % function seq = doGPFA(r, dtMS, nLatents, deleteExistingResults)

            if ~exist('deleteExistingResults', 'var')
                deleteExistingResults = false;
            end
            
            % if sequence data is not yet loaded, do so now
            if isempty(r.sequenceData)
                r.loadSequenceData();
            end
            seqs = r.sequenceData;
            
            % need the training and validation inds
            if isempty(r.inputInfo)
                r.loadInputInfo();
            end
            
            trainInds = r.inputInfo.trainInds;
            validInds = r.inputInfo.validInds;

            xDim = nLatents;

            gpfaOutputFile = r.getGpfaResultsFile(dtMS, nLatents);

            % delete existing file if requested
            if deleteExistingResults
                if exist(gpfaOutputFile, 'file')
                    warning('Deleting existing GPFA results file %s', gpfaOutputFile);
                    delete(gpfaOutputFile)
                end
            end

            % load existing results if they exist
            if exist(gpfaOutputFile, 'file')
                fprintf('Loading previously-run GPFA results (this may take some time)\n');
                resultsOut = r.loadGPFA(dtMS, nLatents);
            else
                resultsOut = cell(numel(seqs), 1);
                binWarned = false;
                prog = ProgressBar(numel(seqs), 'Running GPFA on each dataset');
                for nseq = 1:numel(seqs)
                    prog.update(nseq);
                    seq = seqs{nseq};

                    rebin = dtMS / seq(1).params.dtMS;
                    if rebin ~= dtMS && ~binWarned
                        fprintf(['doGPFA: FYI, sequence data is binned ' ...
                                      'at %i ms, adjusting GPFA binning ' ...
                                      'appropriately, but this is ' ...
                                      'probably dumb\n'], seq(1).dtMS);
                        binWarned = true;
                    end

                    % GPFA is expecting data to be in a field
                    % called 'spikes'
                    [seq.spikes] = seq.y;
                    seq = rmfield(seq, 'y');

                    % GPFA is expecting each trial to be numbered
                    % with a "trialId"
                    tid = num2cell(1:numel(seq), 1, ones(numel(seq), ...
                                                      1));
                    [seq.trialId] = tid{:};
                    tic;
                    % run GPFA with a modified version that allows you to
                    % specify an output directory

                    thisSetGpfaResultsDir = fullfile(r.getGpfaResultsDir(dtMS, nLatents), ...
                                                     sprintf('set%03i', ...
                                                             nseq));
                    results = GPFA.neuralTraj_mod(1, seq, 'trainTrialIdx', trainInds{nseq}, ...
                                             'testTrialIdx', validInds{nseq}, 'saveDir', thisSetGpfaResultsDir,...
                                             'xDim', xDim, 'binWidth', rebin);
                    % postprocess - orthonormalization and "cleanup" (?)
                    [estParams2, seqTrain2, seqTest2] = GPFA.Util.postprocess(results);

                    % trim some of the unnecessary fields of the seqTrains
                    trains = {'seqTrain2','seqTest2'};
                    f2keep = {'trialId','T','y','xsm','xorth'};
                    for nn = 1:numel(trains)
                        tr = eval(trains{nn});
                        fs = fields(tr);
                        f2remove = setdiff(fs, f2keep);
                        for nf = 1:numel(f2remove)
                            tr = rmfield(tr, f2remove{nf});
                        end
                        eval(sprintf('%s = tr;', trains{nn}));
                    end

                    results.seqTrain = seqTrain2;
                    results.seqTest = seqTest2;
                    results.train_inds = trainInds{nseq};
                    results.valid_inds = validInds{nseq};
                    results.estParams = estParams2;
                    toc;

                    resultsOut{nseq} = results;
                end

                % check for necessity of v7.3
                varinfo=whos('resultsOut');
                saveopt='';
                if varinfo.bytes >= 2^31
                    saveopt='-v7.3';
                end

                save(gpfaOutputFile, 'resultsOut', ...
                     saveopt);
                prog.finish();
            end
           

            resultsOut = makecol(resultsOut);
            r.gpfaSequenceData = resultsOut;
        end
        
        function tdSet = addGPFAResultsToTrialData(r)
            % function seqs = addGPFAResultsToSeq(r)
            % returns a sequence that has posterior mean
            % values integrated

            if isempty(r.gpfaSequenceData)
                error(['first load results ' ...
                       'using r.doGPFA( ... )']);
            end

            if isempty(r.trialDataSet)
                r.loadTrialDataFromDatasetCollection();
            end
            tdSet = r.trialDataSet; 
            
            r.loadPosteriorMeans();

            % iterate over datasets
            prog = ProgressBar(r.nDatasets, 'Adding GPFA results to TrialData');
            for iDS = 1:r.nDatasets
                prog.update(iDS);
                gs = r.gpfaSequenceData{iDS};
                
                td = tdSet{iDS};
                %pm = r.posteriorMeans(iDS);
                
                % data must be nTrials x nTime x nChannels tensor
                nGP = size(gs.seqTrain(1).xorth, 1);
                nTime = size(gs.seqTrain(1).xorth, 2);
                gpfa_xorth = nan(td.nTrials, nTime, nGP);
                
                binMs = gs.binWidth;
                preKeep = double(r.params.preKeep);
                postKeep = double(r.params.postKeep);
                gpfa_time = (-preKeep:binMs:postKeep-binMs)';
                assert(numel(gpfa_time) == nTime);

                % training data
                for itr = 1:numel(gs.seqTrain)
                    ntr = gs.seqTrain(itr).trialId;
                    gpfa_xorth(ntr, :, :) = gs.seqTrain(itr).xorth';
                end

                % test data
                for itr = 1:numel(gs.seqTest)
                    ntr = gs.seqTest(itr).trialId;
                    gpfa_xorth(ntr, :, :) = gs.seqTest(itr).xorth';
                end

                td = td.dropAnalogChannelGroup({'gpfa_xorth'}); 
                td = td.addAnalogChannelGroup('gpfa_xorth', genNames('gpfa', nGP), ...
                    gpfa_xorth, gpfa_time, 'timeField', 'gpfa_time', 'isAligned', true);
                
                tdSet{iDS} = td;
            end
            prog.finish();
            
            tdSet = r.setupTrialDataAsLFADS(tdSet);
            r.trialDataSet = tdSet;
            
            function names = genNames(pre, n)
                names = arrayfun(@(i) sprintf('%s%i', pre, i), (1:n)', 'UniformOutput', false);
            end
        end
    end
end
