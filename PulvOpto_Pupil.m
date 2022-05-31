clear

addpath('C:\Users\angel\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Ephys\');
skipRec = 1;

animalCode = '0168';

PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
fileInfo   = dir([PreprocessDir animalCode '_Opto*']); % detect files to load/convert  '_LateralVideo*'


% loop through each recording
for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    %recName = '0168_Opto_010_20180713';
    splitName   = strsplit(recName,'_');
    if datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') < datetime('20180808', 'InputFormat', 'yyyyMMdd'); continue;end

rootPreprocessDir = ['D:\FerretData\' animalCode '\Preprocessed\' recName '\'];
rootAnalysisDir   = ['D:\FerretData\' animalCode '\Analyzed\' recName '\Pupil\'];
rootBehavDatDir   = ['D:\FerretData\' animalCode '\behav\' recName];
if exist(join(rootAnalysisDir),'dir') % skip already analyzed records
    fprintf('Record %s already analyzed \n',recName'); 
    if skipRec == 1; continue; end; end

% load and process ttl data in ephys
load([rootPreprocessDir 'triggerData']);
load(rootBehavDatDir)
load([rootPreprocessDir 'adc_data']);
lfpFs = 1000;
% adc channel names
numADC = size(adc_data,1);
if numADC == 4
    numPul = 3;
    adcChn = {'rightX','rightY','rightD'};
elseif numADC == 7
    numPul = 6;
    adcChn = {'rightX','rightY','rightD','leftX','leftY','leftD'};
end

ttlInd = 1;
rawFs = 30000; % USR DEFNE
trialOnset = find(diff(triggerData(ttlInd,:))==1)./rawFs;

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

condNames = {'Left', 'No Stim'};

if datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') < datetime('20180808', 'InputFormat', 'yyyyMMdd')
    condID = [1 4];
    numMaxEvt = 200;
else
    condID = [1];
    numMaxEvt = 50;
end

numConds = numel(condID);

%% extract snippits of pupil measure

evtDat = nan(numConds,numMaxEvt,numPul,numWinsamp);

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
        for iChn = 1:numPul
            
            evtDat(condCount,iEvt,iChn,:) = adc_data(iChn, evtSampWin); %evtSampWin must >0
                       
        end
        catch
        end
        
    end
    
    condCount = condCount + 1;    

end


evtDat_trialAvg = reshape(nanmedian(evtDat,2),[numConds,numPul,numWinsamp]); % avg across events
% don't use squeeze since first dimension might also be 1, but don't want
% to squeeze that
%% plot
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-100)/2*numPul/3 (screensize(4)-150)/2]);
for iChn = 1:numPul
    subplot(3,numPul/3,iChn)
    hold on

    for iCond = flip(1:numConds)
        toPlot = squeeze(evtDat_trialAvg(iCond,iChn,:));
        %toPlot = smoothts(toPlot,'g',3,0.65);
        plot(toPlot, 'LineWidth', 1)
    end
    if iChn == 1
        legend(condNames{1})
        %legend(condNames{2},condNames{1}) % legend(condNames{1},condNames{4}) % NOTE: match plotting order
    end
    title([adcChn{iChn}])
    xlabel('Time [s]');
    ylabel('Amplitude');
    set(gca,'XTick',linspace(1,numWinsamp,5))
    set(gca,'XTickLabel',linspace(twin(1),twin(2),5))
    axis tight
    %ylim([-30 30])

end

if ~exist(join(rootAnalysisDir),'dir'); mkdir(join(rootAnalysisDir)); end
savefig(fig, [rootAnalysisDir 'Pupil_' num2str(twin(1)) '~' num2str(twin(2)) 'sec.fig'],'compact');
saveas(fig, [rootAnalysisDir 'Pupil_' num2str(twin(1)) '~' num2str(twin(2)) 'sec.png']);
end