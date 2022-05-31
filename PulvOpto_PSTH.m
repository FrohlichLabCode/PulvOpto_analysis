clear

addpath('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Ephys\');
skipRec = 1;

%animalCode = '0180';
animalCodes = {'0187','0188'};
for iAnimal = 1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
switch animalCode
    case '0168'
        PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
        AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
        BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
        regionNames = {'LPl','PPC','VC'};
    case '0173'
        PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed/'];
        AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
        BehavDatDir   = ['E:/FerretData/' animalCode '/behav/'];
        regionNames = {'FC','LPl','PPC','VC'};
    case {'0187','0188'}
        baseDir = ['Z:/Individual/Angel/FerretData/'];
        PreprocessDir = [baseDir animalCode '/Preprocessed/'];
        AnalysisDir   = [baseDir animalCode '/Analyzed/'];
        BehavDatDir   = ['Z:/Ferret Data/' animalCode '/behav/'];
        regionNames = {'PPC','VC'};
        allChn = {[1:16],[17:32]};
        ttlInd = 1; 
end
fileInfo   = dir([PreprocessDir animalCode '_Opto*']); % detect files to load/convert  '_LateralVideo*'


% loop through each recording
for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    %if datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180816', 'InputFormat', 'yyyyMMdd'); continue;end

    %recName = '0168_Opto_010_20180713';

rootPreprocessDir = [PreprocessDir recName '\'];
rootAnalysisDir   = [AnalysisDir recName '\PSTH\'];
rootBehavDatDir   = [BehavDatDir recName];
if exist(join(rootAnalysisDir),'dir') % skip already analyzed records
    fprintf('Record %s already analyzed \n',recName'); 
    if skipRec == 1; continue; end; end


% region info
numRegion = numel(regionNames);

% load and process ttl data in ephys
load([rootPreprocessDir 'triggerData']);
load(rootBehavDatDir)
load([rootPreprocessDir 'lfp\lfpValid']);
load([rootPreprocessDir 'validChn']);
files = dir([rootPreprocessDir 'spikes\spk*.mat']);
totalNumChn = length(files); % get total number of channels before exclusion

if strcmp(splitName{3},'01') == 1 || ismember(splitName{1},{'0187','0188'})
    ttlInd = 1;
else
    ttlInd = 2; %trigger is the 2nd row from 2nd trial on for 0168
end
rawFs = 30000; % USR DEFNE
trialOnset = find(diff(triggerData(ttlInd,:))==1)./rawFs;
numMaxEvt = 200;

% declare info about analysis window and binning of PSTH
twin = [-2 2]; % window to analyze (in s)
binSize = 0.020;    % in seconds

%% preprocess session behav data

% behav data column names
Col_TrialNumber = 1;
Col_TrialType = 2;
Col_Completed = 3; % How long the spout light is on prior to trial initiation

condNames = {'Left','Right','Both','No Stim'};
condID    = unique(session_output_data.BehavData(:,Col_TrialType))'; % 1 or [1 4]
numMaxEvt = size(session_output_data.BehavData,1); % 50 or 200
numConds  = numel(condID);

%%
condCount = 1;
for iCond = condID

    analyzeTheseTrials = find( session_output_data.BehavData(:,Col_TrialType) == iCond );  % USR DEFNE
    
    evtTime = trialOnset(analyzeTheseTrials(1:end)); %(1:end)); %(1:14)
    %evtTime = trialOnset;
    display(['computing PSTH ' recName 'condition:' num2str(iCond)]);
    for iChn = 1:totalNumChn
        load([rootPreprocessDir 'spikes\spk_' num2str(iChn)]);
        spks  = spkTime; % spike times in seconds % OLD: ./1000; % convert spike times (ms) to seconds
        [timePSTH,PSTHrate,psthstats,psthTrial] = is_PSTHstats(evtTime,spks,twin,binSize); % CZ: PSTHrate's 1st timept is time of saccade
        PSTH(condCount,iChn,:)  = PSTHrate;
        % normalise to pre saccade firing rate
        preBins = (timePSTH<-0.5); % 50ms before saccade
        frMean  = mean(PSTHrate(preBins));
        frSTD   = std(PSTHrate(preBins));
        frZ(condCount,iChn,:) = (PSTHrate-frMean)/frSTD; % Spike z score
        spkCell{condCount,iChn} = spks; % save spike times in s for later
    end
    
    numBins = numel(timePSTH);
    
    condCount = condCount + 1;
    
end

%% 
data2Analyze = frZ;
numDivX = 5;

ipanel = 1;

screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-100)/2 (screensize(4)-150)/4*numConds]);
for iCond = 1:numConds

    for iRegion = 1:numel(regionNames)
        
        toPlot = squeeze(data2Analyze(iCond,validChn{iRegion},:));
        
        subplot(numConds,numel(regionNames),ipanel)
        
        imagesc(toPlot)
        title(['Z-score FR PSTH: ' regionNames{iRegion} '; ' condNames{condID(iCond)} ])
        xlabel('Time [s]');
        ylabel('Channel');
        set(gca,'XTick',linspace(1,numBins,numDivX))
        set(gca,'XTickLabel',linspace(twin(1),twin(2),numDivX))
        h = colorbar;
        ylabel(h, 'Z-score FR')
        caxis([-0.5 6])
        axis tight
        
        ipanel = ipanel + 1;
    end
    
end
if ~exist(rootAnalysisDir,'dir'); mkdir(rootAnalysisDir); end
savefig(fig, [rootAnalysisDir 'Z-score FR PSTH.fig'],'compact');
saveas(fig, [rootAnalysisDir 'Z-score FR PSTH.png']);


%% chn avged
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-100)/2 (screensize(4)-150)/4]);

for iRegion = 1:numel(regionNames)
    subplot(1,numel(regionNames),iRegion)
    hold on
    
    if length(condID) == 1
        sliceData = reshape(data2Analyze(iCond,validChn{iRegion},:),[numel(validChn{iRegion}),size(data2Analyze,3)]);
        data2Average = sliceData(any(sliceData,2),:);
        toPlot = squeeze(nanmean(data2Average,1));
        toPlot = smoothts(toPlot,'g',3,0.65);
        plot(toPlot, 'LineWidth', 1.5)
        legend(condNames{condID(1)});
    else
        for iCond = flip(1:numConds)
            % delete all channels with all zeros
            sliceData = reshape(data2Analyze(iCond,validChn{iRegion},:),[numel(validChn{iRegion}),size(data2Analyze,3)]);
            data2Average = sliceData(any(sliceData,2),:);
            toPlot = squeeze(nanmean(data2Average,1));
            toPlot = smoothts(toPlot,'g',3,0.65);
            plot(toPlot, 'LineWidth', 1.5)
        end
        legend(condNames{condID(2)}, condNames{condID(1)});% NOTE: match plotting order
    end
    
    title(['Z-score FR PSTH: ' regionNames{iRegion} ])
    xlabel('Time [s]');
    ylabel('Z-score Firing Rate [Hz]');
    set(gca,'XTick',linspace(1,numBins,numDivX))
    set(gca,'XTickLabel',linspace(twin(1),twin(2),numDivX))
    
    
    axis tight
    ylim([-1.5 10])
end

savefig(fig, [rootAnalysisDir 'Z-score FR PSTH_chn-avg.fig'],'compact');
saveas(fig, [rootAnalysisDir 'Z-score FR PSTH_chn-avg.png']);
end
end
