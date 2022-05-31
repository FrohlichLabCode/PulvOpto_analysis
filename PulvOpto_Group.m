clear
clc

addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));

%% parameteres
% region info
regionNames = {'LPl','PPC','VC'};
numRegion   = numel(regionNames);
regionPairs = {[1,2],[1,3],[2,3]};
numRegionPairs = numel(regionPairs);

% Directory info
GroupAnalysisDir = 'D:\FerretData\0168\GroupAnalysis\';

% time info
twin = [-2 2]; %<<<--- interested time window around event
baseTwin = [-1.5 -0.5]; %<<<--- baseline to subtract from spectrogram 
stimTwin = [0 1]; %1sec stimulation, 2-3 sec ITI
tvec = is_load('D:\FerretData\0168\Analyzed\0168_Opto_108_20180817\FC_validChns_new\LPl-PPC\specAll_Left.mat','tvec');
foi  = is_load('D:\FerretData\0168\Analyzed\0168_Opto_108_20180817\FC_validChns_new\LPl-PPC\specAll_Left.mat','foi');
%tvec = twin(1):0.01:twin(2); %make sure tvec length matches avg
basetvecMask = tvec>= baseTwin(1) & tvec<= baseTwin(2);
stimtvecMask = tvec>= stimTwin(1) & tvec<= stimTwin(2);

% Opto info
optoFreqs = [4.3, 11.7, 38]; %theta, alpha, gamma
numFreqs  = numel(optoFreqs);
freqNames = {'Theta','Alpha','Gamma'};
optoFtwins = {[3,6],[10,14],[30,55]}; %define freq band of interest
optoAmps  = [10, 20, 30, 40, 50];
numAmps   = numel(optoAmps);


%%
%dates with 10-50mW parameterization
sessionDates = {'20180808','20180809','20180814','20180815','20180816','20180817','20180820'};
numDates  = numel(sessionDates);


% Initialize vectors to save
sessionAvg.SpecLPl    = nan(numFreqs, numFreqs, numAmps, numDates);%, numel(sponFiles)); % average spectrum of all sessions
sessionAvg.SpecPPC    = nan(numFreqs, numFreqs, numAmps, numDates);%, numel(sponFiles)); % average spectrum of all sessions
sessionAvg.SpecVC     = nan(numFreqs, numFreqs, numAmps, numDates);%, numel(sponFiles)); % average spectrum of all sessions
sessionAvg.PLVLPl_PPC = nan(numFreqs, numFreqs, numAmps, numDates);
sessionAvg.PLVLPl_VC  = nan(numFreqs, numFreqs, numAmps, numDates);
sessionAvg.PLVPPC_VC  = nan(numFreqs, numFreqs, numAmps, numDates);

sessionAvg.SpecLPl_nostim    = nan(numFreqs, numFreqs, numAmps, 2);%, numel(sponFiles)); % average spectrum of all sessions
sessionAvg.SpecPPC_nostim    = nan(numFreqs, numFreqs, numAmps, 2);%, numel(sponFiles)); % average spectrum of all sessions
sessionAvg.SpecVC_nostim     = nan(numFreqs, numFreqs, numAmps, 2);%, numel(sponFiles)); % average spectrum of all sessions
sessionAvg.PLVLPl_PPC_nostim = nan(numFreqs, numFreqs, numAmps, 2);
sessionAvg.PLVLPl_VC_nostim  = nan(numFreqs, numFreqs, numAmps, 2);
sessionAvg.PLVPPC_VC_nostim  = nan(numFreqs, numFreqs, numAmps, 2);



for iDate = 1:numDates    
    sessionDate = sessionDates{iDate};
    
    if sessionDate == '20180808'
        sessionIDs = {{'033','036','039','042','045'};...
                      {'034','037','040','043','046'};...
                      {'035','038','041','044','047'}}; % theta;alpha;gamma sessions
        condID     = [1];
    elseif sessionDate == '20180809'
        sessionIDs = {{'056','060','048','058','052'};...
                      {'055','062','049','057','053'};...
                      {'054','061','050','059','051'}};
        condID     = [1];
    elseif sessionDate == '20180814'
        sessionIDs = {{'068','063','075','070','073'};...
                      {'066','065','076','071','072'};...
                      {'067','064','077','069','074'}};
        condID     = [1];
    elseif sessionDate == '20180815'
        sessionIDs = {{'092','084','079','081','088'};...
                      {'091','086','078','082','089'};...
                      {'090','085','080','083','087'}};
        condID     = [1];
    elseif sessionDate == '20180816'
        sessionIDs = {{'094','098','103','101','105'};...
                      {'093','096','104','100','107'};...
                      {'095','097','102','099','106'}};  
        condID     = [1];
    elseif sessionDate == '20180817'
        sessionIDs = {{'111','109','120','118','116'};...
                      {'113','108','121','119','115'};...
                      {'112','110','122','117','114'}}; % theta;alpha;gamma sessions
        condID     = [1 4];
    elseif sessionDate == '20180820'
        sessionIDs = {{'133','129','124','128','135'};...
                      {'134','130','123','127','137'};...
                      {'132','131','125','126','136'}}; % theta;alpha;gamma sessions
        condID     = [1 4];
    end
                  
      

    %% load funcCon files from same date into cell: stim and noStim 
    for ioptoFreq = 1:numFreqs
        foiMask(ioptoFreq,:) = foi>=optoFtwins{ioptoFreq}(1) & foi<= optoFtwins{ioptoFreq}(2);
        for iAmp = 1:numAmps
            sessionID = sessionIDs{ioptoFreq}{iAmp};
            for iPair = 1:numRegionPairs
                regionPair = [regionNames{regionPairs{iPair}(1)} '-' regionNames{regionPairs{iPair}(2)}]; %eg. 'LPl-PPC'
                if ~strcmp(sessionID, '133') && ~strcmp(sessionID, '132') % incomplete recordings, skip
                stim{ioptoFreq,iAmp,iPair}  = load(['D:\FerretData\0168\Analyzed\0168_Opto_' sessionID '_' sessionDate '\FC_validChns_new\' regionPair '\funcCon_avg_Left.mat']);
                if numel(condID) == 2
                noStim{ioptoFreq,iAmp,iPair}= load(['D:\FerretData\0168\Analyzed\0168_Opto_' sessionID '_' sessionDate '\FC_validChns_new\' regionPair '\funcCon_avg_NoStim.mat']);
                end
                end
            end
        end
    end


    % put average of each session together
    for ioptoFreq = 1:numFreqs % stimulation frequency
        for iFreq = 1:numFreqs % frequency band of interest
            for iAmp = 1:numAmps                
                avgstimPow = nanmean(nanmean(stim{ioptoFreq,iAmp,1}.avgXSpec(foiMask(iFreq,:),stimtvecMask)));
                avgbasePow = nanmean(nanmean(stim{ioptoFreq,iAmp,1}.avgXSpec(foiMask(iFreq,:),basetvecMask))); % get 1 value
                sessionAvg.SpecLPl(ioptoFreq,iFreq,iAmp,iDate) = (avgstimPow - avgbasePow)/avgbasePow;

                avgstimPow = nanmean(nanmean(stim{ioptoFreq,iAmp,1}.avgYSpec(foiMask(iFreq,:),stimtvecMask)));
                avgbasePow = nanmean(nanmean(stim{ioptoFreq,iAmp,1}.avgYSpec(foiMask(iFreq,:),basetvecMask))); % get 1 value
                sessionAvg.SpecPPC(ioptoFreq,iFreq,iAmp,iDate) = (avgstimPow - avgbasePow)/avgbasePow;

                avgstimPow = nanmean(nanmean(stim{ioptoFreq,iAmp,1}.avgPLV(foiMask(iFreq,:),stimtvecMask)));
                avgbasePow = nanmean(nanmean(stim{ioptoFreq,iAmp,1}.avgPLV(foiMask(iFreq,:),basetvecMask))); % get 1 value
                sessionAvg.PLVLPl_PPC(ioptoFreq,iFreq,iAmp,iDate) = (avgstimPow - avgbasePow)/avgbasePow;

                avgstimPow = nanmean(nanmean(stim{ioptoFreq,iAmp,2}.avgYSpec(foiMask(iFreq,:),stimtvecMask)));
                avgbasePow = nanmean(nanmean(stim{ioptoFreq,iAmp,2}.avgYSpec(foiMask(iFreq,:),basetvecMask))); % get 1 value
                sessionAvg.SpecVC(ioptoFreq,iFreq,iAmp,iDate) = (avgstimPow - avgbasePow)/avgbasePow;

                avgstimPow = nanmean(nanmean(stim{ioptoFreq,iAmp,2}.avgPLV(foiMask(iFreq,:),stimtvecMask)));
                avgbasePow = nanmean(nanmean(stim{ioptoFreq,iAmp,2}.avgPLV(foiMask(iFreq,:),basetvecMask))); % get 1 value
                sessionAvg.PLVLPl_VC(ioptoFreq,iFreq,iAmp,iDate) = (avgstimPow - avgbasePow)/avgbasePow;

                avgstimPow = nanmean(nanmean(stim{ioptoFreq,iAmp,3}.avgPLV(foiMask(iFreq,:),stimtvecMask)));
                avgbasePow = nanmean(nanmean(stim{ioptoFreq,iAmp,3}.avgPLV(foiMask(iFreq,:),basetvecMask))); % get 1 value
                sessionAvg.PLVPPC_VC(ioptoFreq,iFreq,iAmp,iDate) = (avgstimPow - avgbasePow)/avgbasePow;


                % calculate the same for no stim trials
                if numel(condID) == 2
                    if sessionDate == '20180817'; ishamDate = 1;
                    elseif sessionDate == '20180820'; ishamDate = 2; end
                avgstimPow = nanmean(nanmean(noStim{ioptoFreq,iAmp,1}.avgXSpec(foiMask(iFreq,:),stimtvecMask)));
                avgbasePow = nanmean(nanmean(noStim{ioptoFreq,iAmp,1}.avgXSpec(foiMask(iFreq,:),basetvecMask))); % get 1 value
                sessionAvg.SpecLPl_nostim(ioptoFreq,iFreq,iAmp,ishamDate) = (avgstimPow - avgbasePow)/avgbasePow;

                avgstimPow = nanmean(nanmean(noStim{ioptoFreq,iAmp,1}.avgYSpec(foiMask(iFreq,:),stimtvecMask)));
                avgbasePow = nanmean(nanmean(noStim{ioptoFreq,iAmp,1}.avgYSpec(foiMask(iFreq,:),basetvecMask))); % get 1 value
                sessionAvg.SpecPPC_nostim(ioptoFreq,iFreq,iAmp,ishamDate) = (avgstimPow - avgbasePow)/avgbasePow;

                avgstimPow = nanmean(nanmean(noStim{ioptoFreq,iAmp,1}.avgPLV(foiMask(iFreq,:),stimtvecMask)));
                avgbasePow = nanmean(nanmean(noStim{ioptoFreq,iAmp,1}.avgPLV(foiMask(iFreq,:),basetvecMask))); % get 1 value
                sessionAvg.PLVLPl_PPC_nostim(ioptoFreq,iFreq,iAmp,ishamDate) = (avgstimPow - avgbasePow)/avgbasePow;

                avgstimPow = nanmean(nanmean(noStim{ioptoFreq,iAmp,2}.avgYSpec(foiMask(iFreq,:),stimtvecMask)));
                avgbasePow = nanmean(nanmean(noStim{ioptoFreq,iAmp,2}.avgYSpec(foiMask(iFreq,:),basetvecMask))); % get 1 value
                sessionAvg.SpecVC_nostim(ioptoFreq,iFreq,iAmp,ishamDate) = (avgstimPow - avgbasePow)/avgbasePow;

                avgstimPow = nanmean(nanmean(noStim{ioptoFreq,iAmp,2}.avgPLV(foiMask(iFreq,:),stimtvecMask)));
                avgbasePow = nanmean(nanmean(noStim{ioptoFreq,iAmp,2}.avgPLV(foiMask(iFreq,:),basetvecMask))); % get 1 value
                sessionAvg.PLVLPl_VC_nostim(ioptoFreq,iFreq,iAmp,ishamDate) = (avgstimPow - avgbasePow)/avgbasePow;

                avgstimPow = nanmean(nanmean(noStim{ioptoFreq,iAmp,3}.avgPLV(foiMask(iFreq,:),stimtvecMask)));
                avgbasePow = nanmean(nanmean(noStim{ioptoFreq,iAmp,3}.avgPLV(foiMask(iFreq,:),basetvecMask))); % get 1 value
                sessionAvg.PLVPPC_VC_nostim(ioptoFreq,iFreq,iAmp,ishamDate) = (avgstimPow - avgbasePow)/avgbasePow;

                end
            end
        end
    end
end

save([GroupAnalysisDir 'sessionAvg_3Freq_5Amp.mat'],'sessionAvg', '-v7.3');

%%
% plot % power change from baseline


for ioptoFreq = 1:numFreqs
    optoName = freqNames{ioptoFreq};
    
    fig = figure('position', [113    49   370   580]);lw = 2; %x,y,width,height
    xLim = [10,50]; yLim = [-0.5,0.9];
    for iFreq = 1:numFreqs
        subplot(3,3,iFreq);
        shadedErrorBar([10:10:50],squeeze(sessionAvg.SpecLPl_nostim(ioptoFreq,iFreq,:,:))',{@mean,@std},{'k-o','markerfacecolor','k'}); hold on;
        shadedErrorBar([10:10:50],squeeze(sessionAvg.SpecLPl(ioptoFreq,iFreq,:,:))',{@mean,@std},{'c-o','markerfacecolor','c'},0.3); 
        title(['LPl ' freqNames{iFreq}]);
        xlim(xLim);ylim(yLim);

        subplot(3,3,3+iFreq);
        shadedErrorBar([10:10:50],squeeze(sessionAvg.SpecPPC_nostim(ioptoFreq,iFreq,:,:))',{@mean,@std},{'k-o','markerfacecolor','k'}); hold on;
        shadedErrorBar([10:10:50],squeeze(sessionAvg.SpecPPC(ioptoFreq,iFreq,:,:))',{@mean,@std},{'c-o','markerfacecolor','c'},0.3); 
        title(['PPC ' freqNames{iFreq}]);xlim(xLim);ylim(yLim);

        subplot(3,3,6+iFreq);
        shadedErrorBar([10:10:50],squeeze(sessionAvg.SpecVC_nostim(ioptoFreq,iFreq,:,:))',{@mean,@std},{'k-o','markerfacecolor','k'}); hold on;
        shadedErrorBar([10:10:50],squeeze(sessionAvg.SpecVC(ioptoFreq,iFreq,:,:))',{@mean,@std},{'c-o','markerfacecolor','c'},0.3); 
        title(['VC ' freqNames{iFreq}]); xlim(xLim);ylim(yLim); 
        xlabel('Amplitude [mW]');    
    end
    savefig(fig, [GroupAnalysisDir 'sessionAvg_' optoName 'Stim_power_3Freq_5Amp.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'sessionAvg_' optoName 'Stim_power_3Freq_5Amp.png']);

    % plot % PLV change from baseline
    fig = figure('position', [ 113    49   370   580]); %x,y,width,height
    lw = 2;
    for iFreq = 1:numFreqs
        subplot(3,3,iFreq);
        shadedErrorBar([10:10:50],squeeze(sessionAvg.PLVLPl_PPC_nostim(ioptoFreq,iFreq,:,:))',{@mean,@std},{'k-o','markerfacecolor','k'}); hold on;
        shadedErrorBar([10:10:50],squeeze(sessionAvg.PLVLPl_PPC(ioptoFreq,iFreq,:,:))',{@mean,@std},{'c-o','markerfacecolor','c'},0.3); 
        title(['LPl-PPC ' freqNames{iFreq}]);
        xlim(xLim);ylim(yLim); 

        subplot(3,3,3+iFreq);
        shadedErrorBar([10:10:50],squeeze(sessionAvg.PLVLPl_VC_nostim(ioptoFreq,iFreq,:,:))',{@mean,@std},{'k-o','markerfacecolor','k'}); hold on;
        shadedErrorBar([10:10:50],squeeze(sessionAvg.PLVLPl_VC(ioptoFreq,iFreq,:,:))',{@mean,@std},{'c-o','markerfacecolor','c'},0.3); 
        title(['LPl-VC ' freqNames{iFreq}]); xlim(xLim);ylim(yLim); 

        subplot(3,3,6+iFreq);
        shadedErrorBar([10:10:50],squeeze(sessionAvg.PLVPPC_VC_nostim(ioptoFreq,iFreq,:,:))',{@mean,@std},{'k-o','markerfacecolor','k'}); hold on;
        shadedErrorBar([10:10:50],squeeze(sessionAvg.PLVPPC_VC(ioptoFreq,iFreq,:,:))',{@mean,@std},{'c-o','markerfacecolor','c'},0.3); 
        title(['PPC-VC ' freqNames{iFreq}]); xlim(xLim);ylim(yLim);  
        xlabel('Amplitude [mW]');    

    end
    
    savefig(fig, [GroupAnalysisDir 'sessionAvg_' optoName 'Stim_PLV_3Freq_5Amp.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'sessionAvg_' optoName 'Stim_PLV_3Freq_5Amp.png']);
end













%% dates with 20mW stim
sessionIDs = {{'012','013','016','018','019','022','024','028','031'};...
              {'011','017','020','025','027','030'};...
              {'026','029','032'}}; % theta;alpha;gamma sessions, all 20mW


condID = [1,4];

%% load files into cell stim
for ioptoFreq = 1:numFreqs
    foiMask(ioptoFreq,:) = foi>=optoFtwins{ioptoFreq}(1) & foi<= optoFtwins{ioptoFreq}(2);
end

for ioptoFreq = 1:numFreqs
    sessions = sessionIDs{ioptoFreq};
    numSession = numel(sessions);
    % Initialize vectors to save
    sessionAvg.SpecLPl    = nan(numFreqs, numFreqs, numSession);%, numel(sponFiles)); % average spectrum of all sessions
    sessionAvg.SpecPPC    = nan(numFreqs, numFreqs, numSession);%, numel(sponFiles)); % average spectrum of all sessions
    sessionAvg.SpecVC     = nan(numFreqs, numFreqs, numSession);%, numel(sponFiles)); % average spectrum of all sessions
    sessionAvg.PLVLPl_PPC = nan(numFreqs, numFreqs, numSession);
    sessionAvg.PLVLPl_VC  = nan(numFreqs, numFreqs, numSession);
    sessionAvg.PLVPPC_VC  = nan(numFreqs, numFreqs, numSession);

    sessionAvg.SpecLPl_nostim    = nan(numFreqs, numFreqs, numSession);%, numel(sponFiles)); % average spectrum of all sessions
    sessionAvg.SpecPPC_nostim    = nan(numFreqs, numFreqs, numSession);%, numel(sponFiles)); % average spectrum of all sessions
    sessionAvg.SpecVC_nostim     = nan(numFreqs, numFreqs, numSession);%, numel(sponFiles)); % average spectrum of all sessions
    sessionAvg.PLVLPl_PPC_nostim = nan(numFreqs, numFreqs, numSession);
    sessionAvg.PLVLPl_VC_nostim  = nan(numFreqs, numFreqs, numSession);
    sessionAvg.PLVPPC_VC_nostim  = nan(numFreqs, numFreqs, numSession);


    for iSession = 1:numel(sessions)
        sessionID = sessions{iSession};
        for iPair = 1:numRegionPairs
            regionPair = [regionNames{regionPairs{iPair}(1)} '-' regionNames{regionPairs{iPair}(2)}]; %eg. 'LPl-PPC'
            sessionFolder = dir(['D:\FerretData\0168\Analyzed\0168_Opto_' sessionID '*']);
            stim{ioptoFreq,iPair}  = load([sessionFolder.folder '\' sessionFolder.name '\FC_validChns_new\' regionPair '\funcCon_avg_Left.mat']);
            noStim{ioptoFreq,iPair}= load([sessionFolder.folder '\' sessionFolder.name '\FC_validChns_new\' regionPair '\funcCon_avg_NoStim.mat']);
        end
        for iFreq = 1:numFreqs % frequency band of interest
            avgstimPow = nanmean(nanmean(stim{ioptoFreq,1}.avgXSpec(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(stim{ioptoFreq,1}.avgXSpec(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.SpecLPl(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            avgstimPow = nanmean(nanmean(stim{ioptoFreq,1}.avgYSpec(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(stim{ioptoFreq,1}.avgYSpec(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.SpecPPC(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            avgstimPow = nanmean(nanmean(stim{ioptoFreq,1}.avgPLV(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(stim{ioptoFreq,1}.avgPLV(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.PLVLPl_PPC(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            avgstimPow = nanmean(nanmean(stim{ioptoFreq,2}.avgYSpec(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(stim{ioptoFreq,2}.avgYSpec(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.SpecVC(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            avgstimPow = nanmean(nanmean(stim{ioptoFreq,2}.avgPLV(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(stim{ioptoFreq,2}.avgPLV(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.PLVLPl_VC(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            avgstimPow = nanmean(nanmean(stim{ioptoFreq,3}.avgPLV(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(stim{ioptoFreq,3}.avgPLV(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.PLVPPC_VC(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            % calculate the same for no stim trials
            avgstimPow = nanmean(nanmean(noStim{ioptoFreq,1}.avgXSpec(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(noStim{ioptoFreq,1}.avgXSpec(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.SpecLPl_nostim(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            avgstimPow = nanmean(nanmean(noStim{ioptoFreq,1}.avgYSpec(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(noStim{ioptoFreq,1}.avgYSpec(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.SpecPPC_nostim(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            avgstimPow = nanmean(nanmean(noStim{ioptoFreq,1}.avgPLV(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(noStim{ioptoFreq,1}.avgPLV(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.PLVLPl_PPC_nostim(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            avgstimPow = nanmean(nanmean(noStim{ioptoFreq,2}.avgYSpec(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(noStim{ioptoFreq,2}.avgYSpec(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.SpecVC_nostim(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            avgstimPow = nanmean(nanmean(noStim{ioptoFreq,2}.avgPLV(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(noStim{ioptoFreq,2}.avgPLV(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.PLVLPl_VC_nostim(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            avgstimPow = nanmean(nanmean(noStim{ioptoFreq,3}.avgPLV(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(noStim{ioptoFreq,3}.avgPLV(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.PLVPPC_VC_nostim(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;        
        end
    end
    save([GroupAnalysisDir 'sessionAvg_' freqNames{ioptoFreq} '_20mW.mat'],'sessionAvg', '-v7.3');
end



%%
% plot % power change from baseline


for ioptoFreq = 1%1:numFreqs
    optoName = freqNames{ioptoFreq};
    colLabel = {'No stim', optoName};
    load([GroupAnalysisDir 'sessionAvg_' optoName '_20mW.mat']);
    
    fig = figure('position', [20    20   600   900]);lw = 2; %x,y,width,height
    yLim = {[-0.3,0.3],[-0.3,1.2],[-0.7,0.8]};
    for iFreq = 1:numFreqs
        subplot(3,3,numFreqs*iFreq-2);
        mat2plot = [squeeze(sessionAvg.SpecLPl_nostim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.SpecLPl(ioptoFreq,iFreq,:))];       
        H = notBoxPlot(mat2plot,'interval','sem');
        % Color 1st column means gray
        %d=[H.data];
        set([H(1).data],'MarkerSize',6,'markerFaceColor',[1,1,1]*0.5,'markerEdgeColor', 'none')
        set([H(1).sdPtch],'FaceColor',[1,1,1]*0.75,'EdgeColor','none')
        set([H(1).semPtch],'FaceColor',[1,1,1]*0.25,'EdgeColor','none')
        set([H(1).mu],'Color','k')
        set([H(2).sdPtch],'FaceColor',[0,0.7,1],'EdgeColor','none')
        set([H(2).semPtch],'FaceColor',[0,0.7,1]*0.7,'EdgeColor','none')
        set([H(2).mu],'Color',[0,0.7,1]*0.4)       
        % add title and legend for shade colors
        title(['LPl ' freqNames{iFreq}]);
        %legend([H(1).sdPtch H(2).sdPtch], 'no stim',[optoName ' stim']); 
        ylabel('% change from baseline');
        xticklabels(colLabel); ylim(yLim{ioptoFreq});

        subplot(3,3,numFreqs*iFreq-1);
        mat2plot = [squeeze(sessionAvg.SpecPPC_nostim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.SpecPPC(ioptoFreq,iFreq,:))];
        H = notBoxPlot(mat2plot,'interval','sem');
        % Color 1st column means gray
        %d=[H.data];
        set([H(1).data],'MarkerSize',6,'markerFaceColor',[1,1,1]*0.5,'markerEdgeColor', 'none')
        set([H(1).sdPtch],'FaceColor',[1,1,1]*0.75,'EdgeColor','none')
        set([H(1).semPtch],'FaceColor',[1,1,1]*0.25,'EdgeColor','none')
        set([H(1).mu],'Color','k')
        set([H(2).sdPtch],'FaceColor',[0,0.7,1],'EdgeColor','none')
        set([H(2).semPtch],'FaceColor',[0,0.7,1]*0.7,'EdgeColor','none')
        set([H(2).mu],'Color',[0,0.7,1]*0.4)
        % add title and legend for shade colors
        title(['PPC ' freqNames{iFreq}]);
        xticklabels(colLabel); ylim(yLim{ioptoFreq});
        % one-way anova comparing no stim and stim, save p values for
        % correcting multiple comparison
               
        subplot(3,3,numFreqs*iFreq);
        mat2plot = [squeeze(sessionAvg.SpecVC_nostim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.SpecVC(ioptoFreq,iFreq,:))];
        H = notBoxPlot(mat2plot,'interval','sem');
        % Color 1st column means gray
        %d=[H.data];
        set([H(1).data],'MarkerSize',6,'markerFaceColor',[1,1,1]*0.5,'markerEdgeColor', 'none')
        set([H(1).sdPtch],'FaceColor',[1,1,1]*0.75,'EdgeColor','none')
        set([H(1).semPtch],'FaceColor',[1,1,1]*0.25,'EdgeColor','none')
        set([H(1).mu],'Color','k')
        set([H(2).sdPtch],'FaceColor',[0,0.7,1],'EdgeColor','none')
        set([H(2).semPtch],'FaceColor',[0,0.7,1]*0.7,'EdgeColor','none')
        set([H(2).mu],'Color',[0,0.7,1]*0.4)
        % add title and legend for shade colors
        title(['VC ' freqNames{iFreq}]);
        xticklabels(colLabel); ylim(yLim{ioptoFreq});

           
    end
    savefig(fig, [GroupAnalysisDir 'sessionAvg_' optoName 'Stim_power_20mW.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'sessionAvg_' optoName 'Stim_power_20mW.png']);
    
    % one-way anova comparing no stim and stim, save p values for correcting multiple comparison
    for iFreq = 1:numFreqs
        mat2plot = [squeeze(sessionAvg.SpecLPl_nostim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.SpecLPl(ioptoFreq,iFreq,:))];       
        p(numFreqs*iFreq-2) = anova1(mat2plot,colLabel);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
        mat2plot = [squeeze(sessionAvg.SpecPPC_nostim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.SpecPPC(ioptoFreq,iFreq,:))];
        p(numFreqs*iFreq-1) = anova1(mat2plot,colLabel);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
        mat2plot = [squeeze(sessionAvg.SpecVC_nostim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.SpecVC(ioptoFreq,iFreq,:))];
        p(numFreqs*iFreq) = anova1(mat2plot,colLabel);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
    end
    p_corrected = bonf_holm(p,.05);
    save([GroupAnalysisDir 'sessionAvg_' optoName 'Stim_power_20mW_pvalue.mat'],'p','p_corrected');
    clear p p_corrected
    
    
    % plot % PLV change from baseline
    fig = figure('position', [20    20   600   900]);lw = 2; %x,y,width,height
    yLim = {[-0.3,0.4],[-0.3,0.6],[-0.35,0.8]}; % diff range for each stim frequency
    for iFreq = 1:numFreqs
        subplot(3,3,numFreqs*iFreq-2);
        
        mat2plot = [squeeze(sessionAvg.PLVLPl_PPC_nostim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.PLVLPl_PPC(ioptoFreq,iFreq,:))];
        % plot boxplot with data points using notBoxPlot.m
        H = notBoxPlot(mat2plot,'interval','sem');
        % Color 1st column means gray
        %d=[H.data];
        set([H(1).data],'MarkerSize',6,'markerFaceColor',[1,1,1]*0.5,'markerEdgeColor', 'none')
        set([H(1).sdPtch],'FaceColor',[1,1,1]*0.75,'EdgeColor','none')
        set([H(1).semPtch],'FaceColor',[1,1,1]*0.25,'EdgeColor','none')
        set([H(1).mu],'Color','k')
        set([H(2).sdPtch],'FaceColor',[0,0.7,1],'EdgeColor','none')
        set([H(2).semPtch],'FaceColor',[0,0.7,1]*0.7,'EdgeColor','none')
        set([H(2).mu],'Color',[0,0.7,1]*0.4)       
        % add title and legend for shade colors
        title(['LPl-PPC ' freqNames{iFreq}]);
        %legend([H(1).sdPtch H(2).sdPtch], 'no stim',[optoName ' stim']); 
        ylabel('% change from baseline');
        xticklabels(colLabel); ylim(yLim{ioptoFreq});


        subplot(3,3,numFreqs*iFreq-1);
        mat2plot = [squeeze(sessionAvg.PLVLPl_VC_nostim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.PLVLPl_VC(ioptoFreq,iFreq,:))];
        % plot boxplot with data points using notBoxPlot.m
        H = notBoxPlot(mat2plot,'interval','sem');
        % Color 1st column means gray
        %d=[H.data];
        set([H(1).data],'MarkerSize',6,'markerFaceColor',[1,1,1]*0.5,'markerEdgeColor', 'none')
        set([H(1).sdPtch],'FaceColor',[1,1,1]*0.75,'EdgeColor','none')
        set([H(1).semPtch],'FaceColor',[1,1,1]*0.25,'EdgeColor','none')
        set([H(1).mu],'Color','k')
        set([H(2).sdPtch],'FaceColor',[0,0.7,1],'EdgeColor','none')
        set([H(2).semPtch],'FaceColor',[0,0.7,1]*0.7,'EdgeColor','none')
        set([H(2).mu],'Color',[0,0.7,1]*0.4)
        % add title and legend for shade colors
        title(['LPl-VC' freqNames{iFreq}]);
        xticklabels(colLabel); ylim(yLim{ioptoFreq});

        
        subplot(3,3,numFreqs*iFreq);
        mat2plot = [squeeze(sessionAvg.PLVPPC_VC_nostim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.PLVPPC_VC(ioptoFreq,iFreq,:))];
        % plot boxplot with data points using notBoxPlot.m
        H = notBoxPlot(mat2plot,'interval','sem');
        % Color 1st column means gray
        d=[H.data];
        set([H(1).data],'MarkerSize',6,'markerFaceColor',[1,1,1]*0.5,'markerEdgeColor', 'none')
        set([H(1).sdPtch],'FaceColor',[1,1,1]*0.75,'EdgeColor','none')
        set([H(1).semPtch],'FaceColor',[1,1,1]*0.25,'EdgeColor','none')
        set([H(1).mu],'Color','k')
        set([H(2).sdPtch],'FaceColor',[0,0.7,1],'EdgeColor','none')
        set([H(2).semPtch],'FaceColor',[0,0.7,1]*0.7,'EdgeColor','none')
        set([H(2).mu],'Color',[0,0.7,1]*0.4)
        % add title and legend for shade colors
        title(['PPC-VC ' freqNames{iFreq}]);
        xticklabels(colLabel); ylim(yLim{ioptoFreq});

        
    end
    savefig(fig, [GroupAnalysisDir 'sessionAvg_' optoName 'Stim_PLV_20mW.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'sessionAvg_' optoName 'Stim_PLV_20mW.png']);
    
    % one-way anova comparing no stim and stim, save p values for correcting multiple comparison 
    for iFreq = 1:numFreqs    
        mat2plot = [squeeze(sessionAvg.PLVLPl_PPC_nostim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.PLVLPl_PPC(ioptoFreq,iFreq,:))];
        p(numFreqs*iFreq-2) = anova1(mat2plot,colLabel);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
        mat2plot = [squeeze(sessionAvg.PLVLPl_VC_nostim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.PLVLPl_VC(ioptoFreq,iFreq,:))];
        p(numFreqs*iFreq-1) = anova1(mat2plot,colLabel);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
        mat2plot = [squeeze(sessionAvg.PLVPPC_VC_nostim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.PLVPPC_VC(ioptoFreq,iFreq,:))];
        p(numFreqs*iFreq) = anova1(mat2plot,colLabel);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
    end
    p_corrected = bonf_holm(p,.05);
    save([GroupAnalysisDir 'sessionAvg_' optoName 'Stim_PLV_20mW_pvalue.mat'],'p','p_corrected');
    clear p p_corrected
    close all
end
