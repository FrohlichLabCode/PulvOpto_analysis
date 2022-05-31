clear

addpath('C:\Users\angel\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Ephys\');
skipRec = 0;

animalCode = '0173';
switch animalCode
    case '0168'
        PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
        AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
        BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
    case '0173'
        PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed/'];
        AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
        BehavDatDir   = ['E:/FerretData/' animalCode '/behav/'];
end
fileInfo   = dir([PreprocessDir animalCode '_Opto*']); % detect files to load/convert  '_LateralVideo*'


% loop through each recording
for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    %recName = '0168_Opto_010_20180713';
    splitName   = strsplit(recName,'_');
    %if datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180712', 'InputFormat', 'yyyyMMdd'); continue;end

    rootPreprocessDir = [PreprocessDir recName '\'];
    rootAnalysisDir   = [AnalysisDir recName '\ERP\'];
    rootBehavDatDir   = [BehavDatDir recName];
    if exist(join(rootAnalysisDir),'dir') % skip already analyzed records
        fprintf('Record %s already analyzed \n',recName'); 
        if skipRec == 1; continue; end; end

% region info
if animalCode == '0173'
    regionNames = {'FC','LPl','PPC','VC'};
elseif animalCode == '0168'
    regionNames = {'LPl','PPC','VC'};
end
numRegion = numel(regionNames);


% load and process ttl data in ephys
load([rootPreprocessDir 'triggerData']);
load(rootBehavDatDir)
load([rootPreprocessDir 'lfp\lfpValid']);

lfpMat = lfp.reorderedSig;
lfpFs  = lfp.Fs;
allChn = lfp.reorderedChn;
validNumChn = numel([lfp.reorderedChn{:}]);

ttlInd = 1;
rawFs = 30000; % USR DEFNE
trialOnset = find(diff(triggerData(ttlInd,:))==1)./rawFs;
numMaxEvt = 200;

% declare info about analysis window and binning of PSTH
twin = [-0.5 1.5];% window to analyze (in s)
twinSamp = twin*lfpFs;
numWinsamp = numel(twinSamp(1):twinSamp(2)); 
% in seconds

%% preprocess session behav data

% behav data column names
Col_TrialNumber = 1;
Col_TrialType = 2;
Col_Completed = 3; % How long the spout light is on prior to trial initiation

condNames = {'Left','Right','Both','No Stim'};
condID    = unique(session_output_data.BehavData(:,Col_TrialType))'; % 1 or [1 4]
numMaxEvt = size(session_output_data.BehavData,1); % 50 or 200
numConds  = numel(condID);

%% extract snippits of lfp

evtDat = nan(numConds,numMaxEvt,validNumChn,numWinsamp);

condCount = 1;
for iCond = condID
    
    analyzeTheseTrials = find( session_output_data.BehavData(:,Col_TrialType) == iCond );  % USR DEFNE
    if recName(end-3:end) == '0703'
        evtTime = trialOnset(analyzeTheseTrials(1:14)); %%003 session on 0703 (1:14)
    else
        evtTime = trialOnset(analyzeTheseTrials(1:end));
    end
    numEvt = numel(evtTime);

    for iEvt = 1:numEvt
        
        evtSamp = evtTime(iEvt)*lfpFs;
        evtSampWin = round(evtSamp+twinSamp(1):evtSamp+twinSamp(2));
        try
        for iChn = 1:validNumChn
            
            evtDat(condCount,iEvt,iChn,:) = lfpMat(iChn, evtSampWin); %evtSampWin must >0
                       
        end
        catch
        end
        
    end
    
    condCount = condCount + 1;    

end


evtDat_trialAvg = reshape(nanmean(evtDat,2),[numConds,validNumChn,size(evtDat,ndims(evtDat))]); % avg across events
%%
ipanel = 1;
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-100)/2 (screensize(4)-150)/4 * numConds]);
for iCond = 1:numConds
    
    
    for iRegion = 1:numRegion
        
        toPlot = squeeze(evtDat_trialAvg(iCond,allChn{iRegion},:));
        
        subplot(numConds,numel(regionNames),ipanel)
        
        imagesc(toPlot)
        title(['ERP: ' regionNames{iRegion} '; ' condNames{condID(iCond)}])
        xlabel('Time [s]');
        ylabel('Channel');
        set(gca,'XTick',linspace(1,numWinsamp,5))
        set(gca,'XTickLabel',linspace(twin(1),twin(2),5))
        h = colorbar;
        ylabel(h, 'Amplitude [uV]')
        caxis([-20 20])
        axis tight
        
        ipanel = ipanel + 1;
    end
    
end

if ~exist(join(rootAnalysisDir),'dir'); mkdir(join(rootAnalysisDir)); end
savefig(fig, [rootAnalysisDir 'ERP_' num2str(twin(1)) '~' num2str(twin(2)) 'sec.fig'],'compact');
saveas(fig, [rootAnalysisDir 'ERP_' num2str(twin(1)) '~' num2str(twin(2)) 'sec.png']);
%% chn avged
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-100)/2 (screensize(4)-150)/4]);
for iRegion = 1:numel(regionNames)
    subplot(1,numel(regionNames),iRegion)
    hold on
    
    if length(condID) == 1
        toPlot = squeeze(nanmedian(evtDat_trialAvg(1,allChn{iRegion},:),2));
        toPlot = smoothts(toPlot,'g',3,0.65);
        numBins = size(toPlot,2);
        plot(toPlot, 'LineWidth', 1)
        legend(condNames{condID(1)});
    else
        for iCond = flip(1:numConds) % control in the back
        toPlot = squeeze(nanmedian(evtDat_trialAvg(iCond,allChn{iRegion},:),2));
        toPlot = smoothts(toPlot,'g',3,0.65);
        numBins = size(toPlot,2);
        plot(toPlot, 'LineWidth', 1)
        end
        legend(condNames{condID(2)}, condNames{condID(1)});% NOTE: match plotting order
    end

    title(['ERP: ' regionNames{iRegion} ])
    xlabel('Time [s]');
    ylabel('Amplitude [uV]');
    set(gca,'XTick',linspace(1,numWinsamp,5))
    set(gca,'XTickLabel',linspace(twin(1),twin(2),5))

    axis tight
    ylim([-30 30])

end
savefig(fig, [rootAnalysisDir 'ERP_chn-avg_' num2str(twin(1)) '~' num2str(twin(2)) 'sec.fig'],'compact');
saveas(fig, [rootAnalysisDir 'ERP_chn-avg_' num2str(twin(1)) '~' num2str(twin(2)) 'sec.png']);
end