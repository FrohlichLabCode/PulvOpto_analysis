clear

addpath('C:\Users\angel\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Ephys\');
skipRec =1;

% animalCodes = {'0181','0180','0168','0173'};
animalCodes = {'0187','0188'};

for iAnimal = 1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};

    switch animalCode
    case {'0180','0181','0173'}
        PreprocessDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Preprocessed/'];
        AnalysisDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Analyzed/'];
        BehavDatDir   = ['Z:/Ferret Data/' animalCode '/behav/'];
        regionNames = {'PFC','LPl','PPC','VC'};
        allChn = {[1:16],[17:32],[33:48],[49:64]};
        ttlInd = 2;
    case '0168'
        PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
        AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
        BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
        regionNames = {'LPl','PPC','VC'};
        allChn = {[1:16],[17:48],[49:64]};
        ttlInd = 1;
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
    
linORlog = 2;
if linORlog == 1
    numFreqs = 200;
    lowFreq  = 1;
    highFreq = 80;
    foi      = linspace(lowFreq,highFreq,numFreqs); % linear spacing
    fois = [2, 5:5:highFreq];
    tickLabel = string(fois); % generate a string array matches fois {"5","10"...}
elseif linORlog == 2
    numFreqs = 100;
    lowFreq  = 2;
    highFreq = 128;
    foi   = logspace(log10(lowFreq),log10(highFreq),numFreqs); % log spacing
    fois = 2.^(log2(lowFreq):1:log2(highFreq)); %[2 4 8 12 16 32 64 128];
    tickLabel = string(fois);
end

for fi = 1:numel(fois)
    [bi,bb] = sort(abs(foi-fois(fi)));
    tickLoc(fi) = bb(1);
end


% loop through each recording
for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    %recName = '0168_Opto_010_20180713';
    splitName = strsplit(recName,'_');
    if datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180724', 'InputFormat', 'yyyyMMdd'); continue;end

rootPreprocessDir = [PreprocessDir recName '/'];
rootAnalysisDir   = [AnalysisDir recName '/LFP/'];
rootBehavDatDir   = [BehavDatDir recName];
if exist(join(rootAnalysisDir),'dir') % skip already analyzed records
    fprintf('Record %s already analyzed \n',recName'); 
    if skipRec == 1; continue; end; end

% region info
numRegion = numel(regionNames);

% load and process ttl data in ephys
display(['Loading TTL and LFP for rec ' recName]);
load([rootPreprocessDir 'triggerData']);
load(rootBehavDatDir)
load([rootPreprocessDir 'lfp\lfpMat']); % including lfpFs=1000

if splitName{4}(end-3:end) == '0712'
    allChn = {[49:64],[17:48],[1:16]};
end
validNumChn = numel(allChn);
% load([rootPreprocessDir 'lfp\lfpValid']);
% lfpMat = lfp.reorderedSig;
% lfpFs  = lfp.Fs;
% allChn = lfp.reorderedChn;
% validNumChn = numel([lfp.reorderedChn{:}]);

rawFs = 30000; % USR DEFNE
trialOnset = find(diff(triggerData(ttlInd,:))==1)./rawFs;

% % declare info about analysis window and binning of PSTH
% twin = [-0.5 1.5];% window to analyze (in s)
% twinSamp = twin*lfpFs;
% numWinsamp = numel(twinSamp(1):twinSamp(2)); 
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
condCount = 1;
for iCond = condID
    
    analyzeTheseTrials = find( session_output_data.BehavData(:,Col_TrialType) == iCond );  % USR DEFNE
    numTrigger = length(trialOnset);
    if recName(end-3:end) == '0703'
        evtTime = trialOnset(analyzeTheseTrials(1:14)); %%003 session on 0703 (1:14)
    else
        evtTime = trialOnset(analyzeTheseTrials(analyzeTheseTrials<=numTrigger));
    end
    
    numEvt = numel(evtTime);
    xRange = [evtTime(2)-2,evtTime(2)+10]; % in sec
    tRange = round(xRange(1)*lfpFs):round(xRange(end)*lfpFs); %column index
    for iRegion = 1:numRegion
        nChn = numel(allChn{iRegion});
        nCol = 4;
        nRow = nChn/nCol;
        screensize = get( groot, 'Screensize' );
        fig = figure('Position',[10 50 (screensize(3)-100)*nCol/4 (screensize(4)-150)*nRow/8]);
        for iChn = 1:nChn
            subplot(nRow,nCol,iChn)
            plot(tRange/lfpFs, lfpMat(allChn{iRegion}(iChn),tRange));
            vline(evtTime,'k-');
            xlim(xRange);
            title(['chn ' num2str(iChn)]);
            ylim([-200,200]);
            if iChn>(nRow-1)*nCol % last row
                xlabel('Time [s]');
            end
            
            if mod(iChn,nCol)==0 
                set(gca,'yticklabel',{[]})
                ylabel('')
            else
                ylabel('uV') %left most column has unit
            end
        end
        if ~exist(join(rootAnalysisDir),'dir'); mkdir(join(rootAnalysisDir)); end
        savefig(fig, [rootAnalysisDir regionNames{iRegion} '_LFP.fig'],'compact');
        saveas(fig, [rootAnalysisDir regionNames{iRegion} '_LFP.png']);
    
        % Compute spectrogram of each channel
        window   = 0.5*1024;
        noverlap = round(3*window/4);
        specEpoch     = nan(size(lfpMat(allChn{iRegion}(iChn),tRange),1),numFreqs);

        fig = figure('Position',[10 50 (screensize(3)-100)*nCol/4 (screensize(4)-150)*nRow/8]);
        for iChn = 1:nChn
            subplot(nRow,nCol,iChn)
            [s,f,t] = spectrogram(lfpMat(allChn{iRegion}(iChn),tRange),window,noverlap,foi,lfpFs);
            tvec = t+tRange(1)/lfpFs;
            imagesc(tvec,1:numel(foi),pow2db(abs(s).^2)); %put foi (=f) will plot in linear scale
            %caxis([0 1e7]);
            vline(evtTime,'k-');
            title(['chn ' num2str(iChn)]);
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
            %ylim([1,30]); 
            if iChn>(nRow-1)*nCol % last row
                xlabel('Time [s]');
            end
            
            if mod(iChn,nCol)==1
                ylabel('uV') %left most column has unit
            else
                set(gca,'yticklabel',{[]})
                ylabel('')                
            end
            
            if mod(iChn,nCol)==0
                cl = colorbar('eastoutside'); ylabel(cl,'Power [db]','FontSize',12)        
            end
        end
        colormap(jet)
        if ~exist(join(rootAnalysisDir),'dir'); mkdir(join(rootAnalysisDir)); end
        savefig(fig, [rootAnalysisDir regionNames{iRegion} '_spec.fig'],'compact');
        saveas(fig, [rootAnalysisDir regionNames{iRegion} '_spec.png']);
        
    end       
    %condCount = condCount + 1;    

end
close all % close all figure
end
end