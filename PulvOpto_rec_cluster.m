function PulvOpto_rec_cluster(animalCode,irec, fileInfo,folderSuffix, PreprocessDir, AnalysisDir, BehavDatDir, GroupAnalysisDir,...
        cluster, skipRec, linORlog, MedianorPCA)
    
recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    %if cluster == 0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180726', 'InputFormat', 'yyyyMMdd'); continue;end
    sessionID = splitName{3};
    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/FC' folderSuffix '/']; %eg. FC_median
    rootBehavDatDir   = [BehavDatDir recName];

    if ~exist(join(rootAnalysisDir),'dir') 
        mkdir(join(rootAnalysisDir)); fprintf('\nWorking on record %s =============== \n',recName'); end

    % load and process ttl data in ephys
    lfpMat = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpMat'); % don't feed in denoised data with NaN values
    load([rootPreprocessDir 'triggerData']);
    load(rootBehavDatDir)
    

%%
% Define frequencies of interest. Linear spacing for Phase slope index, and
% logarithmic spacing for all other methods.

% Define frequencies of interest. Linear spacing for Phase slope index, and spacing for all other methods.
if linORlog == 1
    numFreqs = 100;
    lowFreq  = 1;
    highFreq = 80;
    %foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
elseif linORlog == 2
    numFreqs = 150;
    lowFreq  = 2;
    highFreq = 128;
    %foi   = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing
end

% region info
switch animalCode
    case {'0173','0180','0181'}
    regionNames = {'PFC','LPl','PPC','VC'};
    regionPairs = {[1,2],[1,3],[2,3],[2,4],[3,4]};
    case {'0168'}
    regionNames = {'LPl','PPC','VC'};
    regionPairs = {[1,2],[1,3],[2,3]};
end
numRegion   = numel(regionNames);


%% select LFP to load
[lfp.validChn, lfp.reorderedChn] = keepChn(recName);
if MedianorPCA == 0
    %lfpMat = lfpMat; % use the raw lfp signal since GC can't deal with NaN
    lfpFs = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpFs');  
    for i = 1:numel(lfp.validChn) 
        regionChn{i} = lfp.validChn{i}; % Pulvinar, PPC, VC
        regionLFP{i} = lfpMat(lfp.validChn{i},:); % reordered channel correspond to reordered lfp
    end
    
elseif MedianorPCA == 1
    lfpMat = lfp.median;
    lfpFs  = lfp.Fs;
    for i = 1:size(lfp.median,1) %lfp.median is an nChannel by nTimepoint array
        regionChn{i} = i; % Pulvinar, PPC, VC
        regionLFP{i} = lfpMat(i,:); % reordered channel correspond to reordered lfp
    end

elseif MedianorPCA == 2
    lfpMat = lfp.PCA;
    lfpFs  = lfp.Fs;
    for i = 1:size(lfp.PCA,1) %lfp.PCA is an nChannel by nTimepoint array
        regionChn{i} = i; % Pulvinar, PPC, VC
        regionLFP{i} = lfpMat(i,:); % reordered channel correspond to reordered lfp
    end
    
elseif MedianorPCA == 3
    lfpFs = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpFs');
    for i = 1:numel(lfp.validChn)
        regionLFP{i} = lfpMat(lfp.validChn{i}(2),:); %pick the first valid channel from each region
        regionChn{i} = i; % Pulvinar, PPC, VC
    end
end

%% parameteres

twin = [-2 2]; %<<<--- interested time window around event
baseTwin = [-1.5 -0.5]; %<<<--- baseline to subtract from spectrogram 
%tvec = twin(1):twin(2); % time vector in sec

ttlInd = 2;
% if strcmp(animalCode, '0181') || (strcmp(sessionID, '01') ~= 1 && strcmp(animalCode, '0173'))
%     ttlInd = 2;
% else
%     ttlInd = 1;
% end
rawFs = 30000;
trialOnset = find(diff(triggerData(ttlInd,:))==1)./rawFs;


%% preprocess session behav data

% behav data column names
Col_TrialNumber = 1;
Col_TrialType = 2;
Col_Completed = 3; % How long the spout light is on prior to trial initiation

condNames = {'Left','Right','Bilateral','NoStim'};
condID    = unique(session_output_data.BehavData(:,Col_TrialType))'; % 1 or [1 4]
numMaxEvt = size(session_output_data.BehavData,1); % 50 or 200
numConds  = numel(condID);

%% get event times for each condition
for iCond = condID
    
    analyzeTheseTrials = find( session_output_data.BehavData(:,Col_TrialType) == iCond );  % USR DEFNE
    if recName(end-3:end) == '0703'
        evtTime = trialOnset(analyzeTheseTrials(1:14)); %%003 session on 0703 (1:14)
    else
        evtTime = trialOnset(analyzeTheseTrials(1:end));
    end
    evtTime(evtTime < abs(twin(1)) | evtTime > size(lfpMat,2)/lfpFs-twin(2)) = []; % exclude events beyond analysis time window
    
    evtTimes{iCond} = evtTime;
    twins{iCond}    = twin;
    baseTwins{iCond}= baseTwin;
end

% only need to save this for the first time running
save([rootPreprocessDir 'eventTimes'], 'evtTimes', 'condNames', 'condID', 'lfpFs');

%% for each region pairs
regionPair_FunConn(skipRec, linORlog, lowFreq, highFreq, numFreqs, lfpFs,sessionID,...
    evtTimes,twin,baseTwin, condNames, condID, regionPairs, regionNames, ...
    regionLFP, regionChn, rootAnalysisDir,GroupAnalysisDir)

% region_spec_by_trial(skipRec, linORlog, lowFreq, highFreq, numFreqs, lfpFs,...
%     evtTimes,twins,baseTwins, condNames, condID, regionNames, ...
%     regionLFP, regionChn, sessionName, GroupAnalysisDir);
end