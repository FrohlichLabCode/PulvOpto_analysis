clear
clc

addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));

%% parameteres
% region info
regionNames = {'LPl','PPC','VC'};
numRegions  = numel(regionNames);
regionPairs = {[1,2],[1,3],[2,3]};
numRegionPairs = numel(regionPairs);
regionPairNames = {'LPl-PPC','LPl-VC','PPC-VC'};
regionPair_Names= {'LPl_PPC','LPl_VC','PPC_VC'};
regionPairNamesGC = {'LPl->PPC','LPl->VC','PPC->VC';'PPC->LPl','VC->LPl','VC->PPC'};

% Directory info
GroupAnalysisDir = 'E:/FerretData/0168/GroupAnalysis/GC/';

% time info
twin = [-2 2]; %<<<--- interested time window around event
baseTwin = [-1.5 -0.5]; %<<<--- baseline to subtract from spectrogram 
stimTwin = [0 1]; %1sec stimulation, 2-3 sec ITI
%tvec = is_load('E:\FerretData\0168\Analyzed\0168_Opto_108_20180817\FC_validChns_new\LPl-PPC\specAll_Left.mat','tvec');
% tvecGC
tvec = is_load('E:\FerretData\0168\Analyzed\0168_Opto_108_20180817\FC_validChns_new\LPl-PPC\GCAll_Left.mat','tvecGC');
foi  = is_load('E:\FerretData\0168\Analyzed\0168_Opto_108_20180817\FC_validChns_new\LPl-PPC\specAll_Left.mat','foi');
%tvec = twin(1):0.01:twin(2); %make sure tvec length matches avg
basetvecMask = tvec>= baseTwin(1) & tvec<= baseTwin(2);
stimtvecMask = tvec>= stimTwin(1) & tvec<= stimTwin(2);
lowFreq = 2; highFreq = 128; linORlog = 2; 
[tickLoc, tickLabel] = getTickLabel(lowFreq, highFreq, 100, linORlog);

% Opto info
optoFreqs = [4.3, 11.7, 38]; %theta, alpha, gamma
numFreqs  = numel(optoFreqs);
freqNames = {'Theta','Alpha','Gamma'};
optoFtwins = {[3,6],[10,14],[30,55]}; %define freq band of interest
optoAmps  = [10, 20, 30, 40, 50];
numAmps   = numel(optoAmps);
numfoi = length(foi);
numtvec = length(tvec);
for ioptoFreq = 1:numFreqs
   foiMask(ioptoFreq,:) = foi>=optoFtwins{ioptoFreq}(1) & foi<= optoFtwins{ioptoFreq}(2);
end
%%
%dates with 10-50mW parameterization
sessionDates = {'20180808','20180809','20180814','20180815','20180816','20180817','20180820'};
numDates  = numel(sessionDates);


% Initialize vectors to save
sessionAll.GCLPl_PPC = nan(2,numFreqs, numAmps, numDates, numfoi, numtvec); % first dimension is whether it's 1->2 or 2->1
sessionAll.GCLPl_VC  = nan(2,numFreqs, numAmps, numDates, numfoi, numtvec);
sessionAll.GCPPC_VC  = nan(2,numFreqs, numAmps, numDates, numfoi, numtvec);
sessionAll.GCLPl_PPC_noStim = nan(2,numFreqs, numAmps, 2, numfoi, numtvec);
sessionAll.GCLPl_VC_noStim  = nan(2,numFreqs, numAmps, 2, numfoi, numtvec);
sessionAll.GCPPC_VC_noStim  = nan(2,numFreqs, numAmps, 2, numfoi, numtvec);


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
                      {'132','131','125','126','136'}}; 
        condID     = [1 4];
    end
                  
      

    %% load funcCon files from same date into cell: stim and noStim 
    for ioptoFreq = 1:numFreqs
        foiMask(ioptoFreq,:) = foi>=optoFtwins{ioptoFreq}(1) & foi<= optoFtwins{ioptoFreq}(2);
        for iAmp = 1:numAmps
            sessionID = sessionIDs{ioptoFreq}{iAmp};
            for iPair = 1:numRegionPairs
                regionPair = [regionNames{regionPairs{iPair}(1)} '-' regionNames{regionPairs{iPair}(2)}]; %eg. 'LPl-PPC'
                if strcmp(sessionID, '133'); % incomplete recordings
                stim{1,ioptoFreq,iAmp,iPair,iDate}  = is_load(['E:\FerretData\0168\Analyzed\0168_Opto_111_20180817' '\FC_validChns_new\' regionPair '\GCAll_Left.mat'],'avgGC_XtoY');
                stim{2,ioptoFreq,iAmp,iPair,iDate}  = is_load(['E:\FerretData\0168\Analyzed\0168_Opto_111_20180817' '\FC_validChns_new\' regionPair '\GCAll_Left.mat'],'avgGC_YtoX');
                elseif strcmp(sessionID, '132');
                stim{1,ioptoFreq,iAmp,iPair,iDate}  = is_load(['E:\FerretData\0168\Analyzed\0168_Opto_112_20180817' '\FC_validChns_new\' regionPair '\GCAll_Left.mat'],'avgGC_XtoY');
                stim{2,ioptoFreq,iAmp,iPair,iDate}  = is_load(['E:\FerretData\0168\Analyzed\0168_Opto_112_20180817' '\FC_validChns_new\' regionPair '\GCAll_Left.mat'],'avgGC_YtoX');
                else
                stim{1,ioptoFreq,iAmp,iPair,iDate}  = is_load(['E:\FerretData\0168\Analyzed\0168_Opto_' sessionID '_' sessionDate '\FC_validChns_new\' regionPair '\GCAll_Left.mat'],'avgGC_XtoY');
                stim{2,ioptoFreq,iAmp,iPair,iDate}  = is_load(['E:\FerretData\0168\Analyzed\0168_Opto_' sessionID '_' sessionDate '\FC_validChns_new\' regionPair '\GCAll_Left.mat'],'avgGC_YtoX');
                end
                if numel(condID) == 2
                    if sessionDate == '20180817'; ishamDate = 1;
                    elseif sessionDate == '20180820'; ishamDate = 2; end
                    if strcmp(sessionID, '133');  % incomplete recordings
                    noStim{1,ioptoFreq,iAmp,iPair,ishamDate}  = is_load(['E:\FerretData\0168\Analyzed\0168_Opto_111_20180817' '\FC_validChns_new\' regionPair '\GCAll_NoStim.mat'],'avgGC_XtoY');
                    noStim{2,ioptoFreq,iAmp,iPair,ishamDate}  = is_load(['E:\FerretData\0168\Analyzed\0168_Opto_111_20180817' '\FC_validChns_new\' regionPair '\GCAll_NoStim.mat'],'avgGC_YtoX');

                    elseif strcmp(sessionID, '132'); 
                    noStim{1,ioptoFreq,iAmp,iPair,ishamDate}  = is_load(['E:\FerretData\0168\Analyzed\0168_Opto_112_20180817' '\FC_validChns_new\' regionPair '\GCAll_NoStim.mat'],'avgGC_XtoY');
                    noStim{2,ioptoFreq,iAmp,iPair,ishamDate}  = is_load(['E:\FerretData\0168\Analyzed\0168_Opto_112_20180817' '\FC_validChns_new\' regionPair '\GCAll_NoStim.mat'],'avgGC_YtoX');
                    else
                    noStim{1,ioptoFreq,iAmp,iPair,ishamDate}  = is_load(['E:\FerretData\0168\Analyzed\0168_Opto_' sessionID '_' sessionDate '\FC_validChns_new\' regionPair '\GCAll_NoStim.mat'],'avgGC_XtoY');
                    noStim{2,ioptoFreq,iAmp,iPair,ishamDate}  = is_load(['E:\FerretData\0168\Analyzed\0168_Opto_' sessionID '_' sessionDate '\FC_validChns_new\' regionPair '\GCAll_NoStim.mat'],'avgGC_YtoX');
                    end
                end
            end
        end
    end
end
save([GroupAnalysisDir 'GC_stim_noStim_3Freq_5Amp.mat'],'stim','noStim', '-v7.3');

%% merge sessions
for iDate = 1:numDates
    for ioptoFreq = 1:numFreqs % stimulation frequency
        for iAmp = 1:numAmps   % stimulation amplitude
            for iRegionPair = 1:numRegionPairs
                regionPair = regionPair_Names{iRegionPair};
                for iDirection = 1:2
                sessionAll.(['GC' regionPair])(iDirection,ioptoFreq,iAmp,iDate,:,:) = stim{iDirection,ioptoFreq,iAmp,iRegionPair,iDate};
                % calculate the same for no stim trials
                if iDate > 5
                    if iDate == 6; ishamDate = 1;
                    elseif iDate == 7; ishamDate = 2; end
                    sessionAll.(['GC' regionPair '_noStim'])(iDirection,ioptoFreq,iAmp,ishamDate,:,:) = noStim{iDirection,ioptoFreq,iAmp,iRegionPair,ishamDate};
                end
                end
            end
        end
    end
end
save([GroupAnalysisDir 'GC_sessionAll_3Freq_5Amp.mat'],'sessionAll', '-v7.3');

% plotting session average spectrogram (4Hz, 20mW)
for ioptoFreq = 1:numFreqs % stimulation frequency
    freqName = freqNames{ioptoFreq};
    for iAmp = 1:numAmps 
        ampName = num2str(optoAmps(iAmp));
        fig = figure('name',['allSessionMedian_GC_' freqName '_Stim_' ampName 'mW'],'position', [10   20   320*numRegions   270*2]);lw = 2; %x,y,width,height
        iPanel = 1;
        for iRegionPair = 1:numRegionPairs
            regionPairName = regionPair_Names{iRegionPair};
            for iDirection = 1:2
                subplot(2,numRegionPairs,iRegionPair + numRegionPairs*(iDirection-1))
                mat2plot = squeeze(nanmedian(real(sessionAll.(['GC' regionPairName])(iDirection,ioptoFreq,iAmp,:,:,:)),4));
                imagesc(tvec,1:numel(foi),mat2plot);
                xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
                ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                %caxis([-2 2]); 
                cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNamesGC{iDirection,iRegionPair}],'FontSize',12)
            end
        end
        fig.Color='white';
        colormap(jet)
        savefig(fig, [GroupAnalysisDir 'allSessionMedian_GC_' freqName '_Stim_' ampName 'mW.fig'],'compact');
        saveas(fig, [GroupAnalysisDir 'allSessionMedian_GC_' freqName '_Stim_' ampName 'mW.png']);        
    end
end
%mean GC is satuated by outliers, use median


% for no stim sessions, take median from all no stim sessions regardless of
% stim parameter
fig = figure('name','allSessionMedian_GC_noStim','position', [10 20 320*numRegions 270*2]);lw = 2; %x,y,width,height

for iRegionPair = 1:numRegionPairs
    regionPairName = regionPair_Names{iRegionPair};
    for iDirection = 1:2
        subplot(2,numRegionPairs,iRegionPair + numRegionPairs*(iDirection-1))
        mat2plot = squeeze(nanmedian(nanmedian(nanmedian(real(sessionAll.(['GC' regionPairName '_noStim'])(iDirection,:,:,:,:,:)),4),3),2));
        imagesc(tvec,1:numel(foi),mat2plot);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([-2 2]); 
        cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionPairNamesGC{iDirection,iRegionPair}],'FontSize',12)
    end
    fig.Color='white';
    colormap(jet)
end
savefig(fig, [GroupAnalysisDir 'allSessionMedian_GC_no_Stim.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'allSessionMedian_GC_no_Stim.png']);
    


%% calculate time resolved average over FOI(4th dimension)
% to get gradient color
numDiv     = 5; % divide into 8 equal (in samples) parts
CT=cbrewer('seq', 'Blues', numDiv+2);
ColorSet = CT(2:end-1,:); % first 2 colors too light

% time resolved baseline normalized power
for ioptoFreq = 1:numFreqs
    optoName = freqNames{ioptoFreq};
    yLim = {[0.5,1.6],[0.7,2.0],[0.5,2.3]};
    fig = figure('position', [113    49   370   580]);lw = 2; %x,y,width,height

    for iFreq = 1:numFreqs
        subplot(3,3,numFreqs*iFreq-2);
        stimPow = squeeze(nanmean(sessionAll.SpecLPl(ioptoFreq,:,:,foiMask(iFreq,:),:),4)); %numAmps x numDates x t
        shamPow = squeeze(nanmean(nanmean(sessionAll.SpecLPl_noStim(ioptoFreq,:,:,foiMask(iFreq,:),:),4),2));        
        plot(tvec,squeeze(nanmean(shamPow,1)),'k'); hold on;  % combine different amp
        set(gca, 'ColorOrder', ColorSet);
        plot(tvec,squeeze(nanmean(stimPow,2))); % numAmps x t
        ylim(yLim{ioptoFreq}); ylabel([freqNames{iFreq} ' normed power']);
        if iFreq == 1; title(regionNames{1}); %legend('no stim','10mW','20mW','30mW','40mW','50mW'); 
        end
        
        subplot(3,3,numFreqs*iFreq-1);
        stimPow = squeeze(nanmean(sessionAll.SpecPPC(ioptoFreq,:,:,foiMask(iFreq,:),:),4)); %numAmps x numDates x t
        shamPow = squeeze(nanmean(nanmean(sessionAll.SpecPPC_noStim(ioptoFreq,:,:,foiMask(iFreq,:),:),4),2));        
        plot(tvec,squeeze(nanmean(shamPow,1)),'k'); hold on; % combine different amp
        set(gca, 'ColorOrder', ColorSet);
        plot(tvec,squeeze(nanmean(stimPow,2))); % numAmps x t
        ylim(yLim{ioptoFreq});
        if iFreq == 1; title(regionNames{2});end
        
        subplot(3,3,numFreqs*iFreq);
        stimPow = squeeze(nanmean(sessionAll.SpecVC(ioptoFreq,:,:,foiMask(iFreq,:),:),4)); %numAmps x numDates x t
        shamPow = squeeze(nanmean(nanmean(sessionAll.SpecVC_noStim(ioptoFreq,:,:,foiMask(iFreq,:),:),4),2));        
        plot(tvec,squeeze(nanmean(shamPow,1)),'k'); hold on; % combine different amp
        set(gca, 'ColorOrder', ColorSet);
        plot(tvec,squeeze(nanmean(stimPow,2))); % numAmps x t
        ylim(yLim{ioptoFreq});
        if iFreq == 1; title(regionNames{3});end
        if iFreq == 3; xlabel('Time [s]');end
    end
    savefig(fig, [GroupAnalysisDir 'timeResolved_' optoName 'Stim_power_3Freq_5Amp.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'timeResolved_' optoName 'Stim_power_3Freq_5Amp.png']); 
end

% time resolved baseline normalized PLV
for ioptoFreq = 1:numFreqs
    optoName = freqNames{ioptoFreq};
    yLim = {[-1.5,3],[-2,4],[-2,4]};
    fig = figure('position', [113    49   370   580]);lw = 2; %x,y,width,height
    

    for iFreq = 1:numFreqs
        subplot(3,3,numFreqs*iFreq-2);
        stimPow = squeeze(nanmean(sessionAll.PLVLPl_PPC(ioptoFreq,:,:,foiMask(iFreq,:),:),4)); %numAmps x numDates x t
        shamPow = squeeze(nanmean(nanmean(sessionAll.PLVLPl_PPC_noStim(ioptoFreq,:,:,foiMask(iFreq,:),:),4),2));
        plot(tvec,squeeze(nanmean(shamPow,1)),'k'); hold on;  % combine different amp
        set(gca, 'ColorOrder', ColorSet);
        plot(tvec,squeeze(nanmean(stimPow,2))); % numAmps x t
        ylim(yLim{ioptoFreq}); 
        ylabel([freqNames{iFreq} ' PLV normed']);
        if iFreq == 1; title(regionPairNames{1}); %legend('no stim','10mW','20mW','30mW','40mW','50mW'); 
        end
        
        subplot(3,3,numFreqs*iFreq-1);
        stimPow = squeeze(nanmean(sessionAll.PLVLPl_VC(ioptoFreq,:,:,foiMask(iFreq,:),:),4)); %numAmps x numDates x t
        shamPow = squeeze(nanmean(nanmean(sessionAll.PLVLPl_VC_noStim(ioptoFreq,:,:,foiMask(iFreq,:),:),4),2));        
        plot(tvec,squeeze(nanmean(shamPow,1)),'k'); hold on; % combine different amp
        set(gca, 'ColorOrder', ColorSet);
        plot(tvec,squeeze(nanmean(stimPow,2))); % numAmps x t
        ylim(yLim{ioptoFreq});
        if iFreq == 1; title(regionPairNames{2});end
        
        subplot(3,3,numFreqs*iFreq);
        stimPow = squeeze(nanmean(sessionAll.PLVPPC_VC(ioptoFreq,:,:,foiMask(iFreq,:),:),4)); %numAmps x numDates x t
        shamPow = squeeze(nanmean(nanmean(sessionAll.PLVPPC_VC_noStim(ioptoFreq,:,:,foiMask(iFreq,:),:),4),2));        
        plot(tvec,squeeze(nanmean(shamPow,1)),'k'); hold on; % combine different amp
        set(gca, 'ColorOrder', ColorSet);
        plot(tvec,squeeze(nanmean(stimPow,2))); % numAmps x t
        ylim(yLim{ioptoFreq});
        if iFreq == 1; title(regionPairNames{3});end
        if iFreq == 3; xlabel('Time [s]');end
    end
    savefig(fig, [GroupAnalysisDir 'timeResolved_' optoName 'Stim_PLV_3Freq_5Amp.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'timeResolved_' optoName 'Stim_PLV_3Freq_5Amp.png']);         
end

%% time resolved baseline normalized PLV + sem
for ioptoFreq = 1:numFreqs
    optoName = freqNames{ioptoFreq};
    yLim = {[-3,5],[-3,5],[-3,5]};
    fig = figure('name',[optoName 'Stim_tPLV+sem_3Freq_5Amp'], 'position', [10 20 650 900]);lw = 4; %x,y,width,height [10 20 370 580]
    
for iRegionPair = 1:numRegionPairs
    regionPair  = regionPair_Names{iRegionPair};
    for iFreq = 1:numFreqs
        subplot(3,3,numFreqs*iFreq-3+iRegionPair);
        stimPow = squeeze(nanmean(sessionAll.(['PLV' regionPair])(ioptoFreq,:,:,foiMask(iFreq,:),:),4)); %numAmps x numDates x t
        shamPow = squeeze(nanmean(nanmean(sessionAll.(['PLV' regionPair '_noStim'])(ioptoFreq,:,:,foiMask(iFreq,:),:),4),2));
        mean = squeeze(nanmean(shamPow,1));
        sem  = squeeze(nanstd(shamPow,[],1)/sqrt(size(shamPow,1)));
        shadedErrorBar(tvec,mean,sem,{'k','LineWidth',lw}); hold on;  % combine different amp
        set(gca, 'ColorOrder', ColorSet);
        for iAmp = 1:numAmps
            mean = squeeze(nanmean(stimPow(iAmp,:,:),2));
            sem  = squeeze(nanstd(stimPow(iAmp,:,:),[],2));
            shadedErrorBar(tvec,mean,sem,{'color',ColorSet(iAmp,:),'LineWidth',lw}, 0.05);hold on;
        end
        ylim(yLim{ioptoFreq}); xlim([-1,2]);
        % get rid of some labels b/c this figures can't be processed in Adobe
        if iRegionPair ==1; ylabel([freqNames{iFreq} ' PLV normed']);
        else yticklabels([]); end
        if iFreq == 1; title(regionPairNames{iRegionPair});xticklabels([]); %legend('no stim','10mW','20mW','30mW','40mW','50mW'); 
        elseif iFreq == 3; xlabel('Time [s]');
        else xticklabels([]); end
    end
    fig.Color = 'white';
    savefig(fig, [GroupAnalysisDir 'SessionMean_tPLV+sem_' optoName 'Stim_3Freq_5Amp.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'SessionMean_tPLV+sem_' optoName 'Stim_3Freq_5Amp.png']);         
end
end






%% get average of FOI and time for each session -- also for stats
for iFreq = 1:numFreqs
    sessionAvg.SpecLPl(:,iFreq,:,:) = squeeze(nanmean(nanmean(sessionAll.SpecLPl(:,:,:,foiMask(iFreq,:),stimtvecMask),4),5)); %numFreqs x numAmps x numDates
    sessionAvg.SpecPPC(:,iFreq,:,:) = squeeze(nanmean(nanmean(sessionAll.SpecPPC(:,:,:,foiMask(iFreq,:),stimtvecMask),4),5)); %numFreqs x numAmps x numDates
    sessionAvg.SpecVC(:,iFreq,:,:)  = squeeze(nanmean(nanmean(sessionAll.SpecVC(:,:,:,foiMask(iFreq,:),stimtvecMask),4),5)); %numFreqs x numAmps x numDates
    sessionAvg.SpecLPl_noStim(:,iFreq,:,:) = squeeze(nanmean(nanmean(sessionAll.SpecLPl_noStim(:,:,:,foiMask(iFreq,:),stimtvecMask),4),5));        
    sessionAvg.SpecPPC_noStim(:,iFreq,:,:) = squeeze(nanmean(nanmean(sessionAll.SpecPPC_noStim(:,:,:,foiMask(iFreq,:),stimtvecMask),4),5));        
    sessionAvg.SpecVC_noStim(:,iFreq,:,:)  = squeeze(nanmean(nanmean(sessionAll.SpecVC_noStim(:,:,:,foiMask(iFreq,:),stimtvecMask),4),5));        
    sessionAvg.PLVLPl_PPC(:,iFreq,:,:) = squeeze(nanmean(nanmean(sessionAll.PLVLPl_PPC(:,:,:,foiMask(iFreq,:),stimtvecMask),4),5)); %numFreqs x numAmps x numDates
    sessionAvg.PLVLPl_VC(:,iFreq,:,:)  = squeeze(nanmean(nanmean(sessionAll.PLVLPl_VC(:,:,:,foiMask(iFreq,:),stimtvecMask),4),5)); %numFreqs x numAmps x numDates
    sessionAvg.PLVPPC_VC(:,iFreq,:,:)  = squeeze(nanmean(nanmean(sessionAll.PLVPPC_VC(:,:,:,foiMask(iFreq,:),stimtvecMask),4),5)); %numFreqs x numAmps x numDates
    sessionAvg.PLVLPl_PPC_noStim(:,iFreq,:,:) = squeeze(nanmean(nanmean(sessionAll.PLVLPl_PPC_noStim(:,:,:,foiMask(iFreq,:),stimtvecMask),4),5)); %numFreqs x numAmps x numDates
    sessionAvg.PLVLPl_VC_noStim(:,iFreq,:,:)  = squeeze(nanmean(nanmean(sessionAll.PLVLPl_VC_noStim(:,:,:,foiMask(iFreq,:),stimtvecMask),4),5)); %numFreqs x numAmps x numDates
    sessionAvg.PLVPPC_VC_noStim(:,iFreq,:,:)  = squeeze(nanmean(nanmean(sessionAll.PLVPPC_VC_noStim(:,:,:,foiMask(iFreq,:),stimtvecMask),4),5)); %numFreqs x numAmps x numDates
end
save([GroupAnalysisDir 'sessionAvg_5Amp.mat'],'sessionAvg', '-v7.3');

% % plot avg power over diff amplitudes (too few control trials, can't see sig diff)
% for ioptoFreq = 1:numFreqs
%     optoName = freqNames{ioptoFreq};
%     xLim = [10,50];
%     yLim = {[0.5,2.5],[0.6,3],[0,2.3]};
%     fig = figure('position', [113    49   370   580]);lw = 2; %x,y,width,height
% 
%     for iFreq = 1:numFreqs
%         subplot(3,3,numFreqs*iFreq-2);
%         shadedErrorBar([10:10:50],squeeze(sessionAvg.SpecLPl_noStim(ioptoFreq,iFreq,:,:))',{@mean,@std},{'k-o','markerfacecolor','k'}); hold on;
%         shadedErrorBar([10:10:50],squeeze(sessionAvg.SpecLPl(ioptoFreq,iFreq,:,:))',{@mean,@std},{'-o','markerfacecolor',ColorSet(end,:),'color',ColorSet(end,:)},0.3); 
%         xlim(xLim); %ylim(yLim{ioptoFreq}); ylabel([freqNames{iFreq} ' normed power']);
%         if iFreq == 1; title(regionNames{1}); %legend('no stim','10mW','20mW','30mW','40mW','50mW'); 
%         end
%         
%         subplot(3,3,numFreqs*iFreq-1);
%         shadedErrorBar([10:10:50],squeeze(sessionAvg.SpecPPC_noStim(ioptoFreq,iFreq,:,:))',{@mean,@std},{'k-o','markerfacecolor','k'}); hold on;
%         shadedErrorBar([10:10:50],squeeze(sessionAvg.SpecPPC(ioptoFreq,iFreq,:,:))',{@mean,@std},{'-o','markerfacecolor',ColorSet(end,:),'color',ColorSet(end,:)},0.3);         
%         %ylim(yLim{ioptoFreq});
%         if iFreq == 1; title(regionNames{2});end
%         
%         subplot(3,3,numFreqs*iFreq);
%         shadedErrorBar([10:10:50],squeeze(sessionAvg.SpecVC_noStim(ioptoFreq,iFreq,:,:))',{@mean,@std},{'k-o','markerfacecolor','k'}); hold on;
%         shadedErrorBar([10:10:50],squeeze(sessionAvg.SpecVC(ioptoFreq,iFreq,:,:))',{@mean,@std},{'-o','markerfacecolor',ColorSet(end,:),'color',ColorSet(end,:)},0.3); 
%         %ylim(yLim{ioptoFreq});
%         if iFreq == 1; title(regionNames{3});end
%         if iFreq == 3; xlabel('Time [s]');end
%     end
%     savefig(fig, [GroupAnalysisDir 'sessionAvg_' optoName 'Stim_power_5Amp.fig'],'compact');
%     saveas(fig, [GroupAnalysisDir 'sessionAvg_' optoName 'Stim_power_5Amp.png']); 
% end

%% box plot and stats on 0v10v20v30v40v50 mA
colLabel = {'No stim', '10', '20','30','40','50mW'};
groupVec = [zeros(1,length(noStimAll)),reshape(repmat([1;2;3;4;5],1,numDates)',1,[])]; %eg. 000111122223333

%------------- for normalized power
for ioptoFreq = 1:numFreqs
    optoName = freqNames{ioptoFreq};

    fig = figure('name', 'boxPlot_power_5Amp', 'position', [113    49   370   580]);lw = 2; %x,y,width,height
    yLim = {[0.5,2.5],[0.5,2.5],[0.5,2.5]};
    for iFreq = 1:numFreqs
        subplot(3,3,numFreqs*iFreq-2);
        noStimAll = reshape(sessionAvg.SpecLPl_noStim(ioptoFreq,iFreq,:,:),1,[]); % collapse all amplitude 
        mat2plot{ioptoFreq,iFreq,1} = [noStimAll, reshape(squeeze(sessionAvg.SpecLPl(ioptoFreq,iFreq,:,:))',1,[])]; % combine 10mA(x7trials), 20mA ... into 1 row
        boxplot(mat2plot{ioptoFreq,iFreq,1},groupVec);
        sub_color_myBars(CT, numDiv+1);
        ylim(yLim{ioptoFreq}); 
        %ylabel([freqNames{iFreq} ' power normed']); % add ylabel distorts figure shape
        if iFreq == 1; title(regionNames{1}); %legend('no stim','10mW','20mW','30mW','40mW','50mW'); 
        end
        
        subplot(3,3,numFreqs*iFreq-1);
        noStimAll = reshape(sessionAvg.SpecPPC_noStim(ioptoFreq,iFreq,:,:),1,[]); % collapse all amplitude 
        mat2plot{ioptoFreq,iFreq,2} = [noStimAll, reshape(squeeze(sessionAvg.SpecPPC(ioptoFreq,iFreq,:,:))',1,[])]; % combine 10mA(x7trials), 20mA ... into 1 row
        boxplot(mat2plot{ioptoFreq,iFreq,2},groupVec);
        sub_color_myBars(CT, numDiv+1);
        ylim(yLim{ioptoFreq}); 
        if iFreq == 1; title(regionNames{2});end
        
        subplot(3,3,numFreqs*iFreq);
        noStimAll = reshape(sessionAvg.SpecVC_noStim(ioptoFreq,iFreq,:,:),1,[]); % collapse all amplitude 
        mat2plot{ioptoFreq,iFreq,3} = [noStimAll, reshape(squeeze(sessionAvg.SpecVC(ioptoFreq,iFreq,:,:))',1,[])]; % combine 10mA(x7trials), 20mA ... into 1 row
        boxplot(mat2plot{ioptoFreq,iFreq,3},groupVec);
        sub_color_myBars(CT, numDiv+1);
        ylim(yLim{ioptoFreq}); 
        if iFreq == 1; title(regionNames{3});end
        %if iFreq == 3; xlabel('Amplitude [mW]');end % add xlabel distorts figure shape
    end
    savefig(fig, [GroupAnalysisDir 'boxPlot_' optoName 'Stim_power_5Amp.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'boxPlot_' optoName 'Stim_power_5Amp.png']); 

end
save([GroupAnalysisDir 'boxPlot_power_5Amp.mat'],'mat2plot', '-v7.3');

for ioptoFreq = 1:numFreqs
    % one-way anova comparing no stim and stim, save p values for correcting multiple comparison
    %p = 3 optoFreq x (3 region x 3 freqs) x (1 main p + 5 pairwise p)   = 3x9x6  
    for iFreq = 1:numFreqs
        [p(ioptoFreq,numFreqs*iFreq-2,1),~,stats] = anova1(mat2plot{ioptoFreq,iFreq,1},groupVec);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
        [c,~,~,gnames] = multcompare(stats);
        p(ioptoFreq,numFreqs*iFreq-2,2:6) = c(1:5,6)'; %7th column is p value, 1-5 row is no stim vs. 10-50mA
        p_corrected(ioptoFreq,numFreqs*iFreq-2,:) = bonf_holm(p(ioptoFreq,numFreqs*iFreq-2,:),.05); % only correct for pair-wise comparison
        
        [p(ioptoFreq,numFreqs*iFreq-1,1),~,stats] = anova1(mat2plot{ioptoFreq,iFreq,2},groupVec);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
        [c,~,~,gnames] = multcompare(stats);
        p(ioptoFreq,numFreqs*iFreq-1,2:6) = c(1:5,6)'; %7th column is p value, 1-5 row is no stim vs. 10-50mA
        p_corrected(ioptoFreq,numFreqs*iFreq-1,:) = bonf_holm(p(ioptoFreq,numFreqs*iFreq-1,:),.05); % only correct for pair-wise comparison

        [p(ioptoFreq,numFreqs*iFreq,1),~,stats] = anova1(mat2plot{ioptoFreq,iFreq,3},groupVec);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
        [c,~,~,gnames] = multcompare(stats);
        p(ioptoFreq,numFreqs*iFreq,2:6) = c(1:5,6)'; %7th column is p value, 1-5 row is no stim vs. 10-50mA
        p_corrected(ioptoFreq,numFreqs*iFreq,:) = bonf_holm(p(ioptoFreq,numFreqs*iFreq,:),.05); % only correct for pair-wise comparison
    
    end
end
save([GroupAnalysisDir 'boxPlot_power_5Amp_pvalue.mat'],'p','p_corrected');
clear p p_corrected

%------------- for normalized PLV
for ioptoFreq = 1:numFreqs
    optoName = freqNames{ioptoFreq};

    fig = figure('position', [113    49   370   580]);lw = 2; %x,y,width,height
    yLim = {[-2,4],[-2,5],[-3,3.3]};
    for iFreq = 1:numFreqs
        subplot(3,3,numFreqs*iFreq-2);
        noStimAll = reshape(sessionAvg.PLVLPl_PPC_noStim(ioptoFreq,iFreq,:,:),1,[]); % collapse all amplitude 
        mat2plot{ioptoFreq,iFreq,1} = [noStimAll, reshape(squeeze(sessionAvg.PLVLPl_PPC(ioptoFreq,iFreq,:,:))',1,[])]; % combine 10mA(x7trials), 20mA ... into 1 row
        boxplot(mat2plot{ioptoFreq,iFreq,1},groupVec);
        sub_color_myBars(CT, numDiv+1);
        %ylim(yLim{ioptoFreq}); 
        %ylabel([freqNames{iFreq} ' power normed']); % add ylabel distorts figure shape
        if iFreq == 1; title(regionPairNames{1}); %legend('no stim','10mW','20mW','30mW','40mW','50mW'); 
        end
        
        subplot(3,3,numFreqs*iFreq-1);
        noStimAll = reshape(sessionAvg.PLVLPl_VC_noStim(ioptoFreq,iFreq,:,:),1,[]); % collapse all amplitude 
        mat2plot{ioptoFreq,iFreq,2} = [noStimAll, reshape(squeeze(sessionAvg.PLVLPl_VC(ioptoFreq,iFreq,:,:))',1,[])]; % combine 10mA(x7trials), 20mA ... into 1 row
        boxplot(mat2plot{ioptoFreq,iFreq,2},groupVec);
        sub_color_myBars(CT, numDiv+1);
        %ylim(yLim{ioptoFreq}); 
        if iFreq == 1; title(regionPairNames{2});end
        
        subplot(3,3,numFreqs*iFreq);
        noStimAll = reshape(sessionAvg.PLVPPC_VC_noStim(ioptoFreq,iFreq,:,:),1,[]); % collapse all amplitude 
        mat2plot{ioptoFreq,iFreq,3} = [noStimAll, reshape(squeeze(sessionAvg.PLVPPC_VC(ioptoFreq,iFreq,:,:))',1,[])]; % combine 10mA(x7trials), 20mA ... into 1 row
        boxplot(mat2plot{ioptoFreq,iFreq,3},groupVec);
        sub_color_myBars(CT, numDiv+1);
        %ylim(yLim{ioptoFreq}); 
        if iFreq == 1; title(regionPairNames{3});end
        %if iFreq == 3; xlabel('Amplitude [mW]');end % add xlabel distorts figure shape
    end
    savefig(fig, [GroupAnalysisDir 'boxPlot_' optoName 'Stim_PLV_5Amp.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'boxPlot_' optoName 'Stim_PLV_5Amp.png']); 
end
save([GroupAnalysisDir 'boxPlot_PLV_5Amp.mat'],'mat2plot', '-v7.3');

for ioptoFreq = 1:numFreqs
    % one-way anova comparing no stim and stim, save p values for correcting multiple comparison
    %p = 3 optoFreq x (3 region x 3 freqs) x (1 main p + 5 pairwise p)   = 3x9x6  
    for iFreq = 1:numFreqs
        [p(ioptoFreq,numFreqs*iFreq-2,1),~,stats] = anova1(mat2plot{ioptoFreq,iFreq,1},groupVec);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
        [c,~,~,gnames] = multcompare(stats);
        p(ioptoFreq,numFreqs*iFreq-2,2:6) = c(1:5,6)'; %7th column is p value, 1-5 row is no stim vs. 10-50mA
        p_corrected(ioptoFreq,numFreqs*iFreq-2,:) = bonf_holm(p(ioptoFreq,numFreqs*iFreq-2,:),.05); % only correct for pair-wise comparison
        
        [p(ioptoFreq,numFreqs*iFreq-1,1),~,stats] = anova1(mat2plot{ioptoFreq,iFreq,2},groupVec);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
        [c,~,~,gnames] = multcompare(stats);
        p(ioptoFreq,numFreqs*iFreq-1,2:6) = c(1:5,6)'; %7th column is p value, 1-5 row is no stim vs. 10-50mA
        p_corrected(ioptoFreq,numFreqs*iFreq-1,:) = bonf_holm(p(ioptoFreq,numFreqs*iFreq-1,:),.05); % only correct for pair-wise comparison

        [p(ioptoFreq,numFreqs*iFreq,1),~,stats] = anova1(mat2plot{ioptoFreq,iFreq,3},groupVec);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
        [c,~,~,gnames] = multcompare(stats);
        p(ioptoFreq,numFreqs*iFreq,2:6) = c(1:5,6)'; %7th column is p value, 1-5 row is no stim vs. 10-50mA
        p_corrected(ioptoFreq,numFreqs*iFreq,:) = bonf_holm(p(ioptoFreq,numFreqs*iFreq,:),.05); % only correct for pair-wise comparison
    
    end
end
save([GroupAnalysisDir 'boxPlot_PLV_5Amp_pvalue.mat'],'p','p_corrected');
clear p p_corrected









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

    sessionAvg.SpecLPl_noStim    = nan(numFreqs, numFreqs, numSession);%, numel(sponFiles)); % average spectrum of all sessions
    sessionAvg.SpecPPC_noStim    = nan(numFreqs, numFreqs, numSession);%, numel(sponFiles)); % average spectrum of all sessions
    sessionAvg.SpecVC_noStim     = nan(numFreqs, numFreqs, numSession);%, numel(sponFiles)); % average spectrum of all sessions
    sessionAvg.PLVLPl_PPC_noStim = nan(numFreqs, numFreqs, numSession);
    sessionAvg.PLVLPl_VC_noStim  = nan(numFreqs, numFreqs, numSession);
    sessionAvg.PLVPPC_VC_noStim  = nan(numFreqs, numFreqs, numSession);


    for iSession = 1:numel(sessions)
        sessionID = sessions{iSession};
        for iPair = 1:numRegionPairs
            regionPair = [regionNames{regionPairs{iPair}(1)} '-' regionNames{regionPairs{iPair}(2)}]; %eg. 'LPl-PPC'
            sessionFolder = dir(['E:\FerretData\0168\Analyzed\0168_Opto_' sessionID '*']);
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
            sessionAvg.SpecLPl_noStim(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            avgstimPow = nanmean(nanmean(noStim{ioptoFreq,1}.avgYSpec(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(noStim{ioptoFreq,1}.avgYSpec(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.SpecPPC_noStim(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            avgstimPow = nanmean(nanmean(noStim{ioptoFreq,1}.avgPLV(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(noStim{ioptoFreq,1}.avgPLV(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.PLVLPl_PPC_noStim(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            avgstimPow = nanmean(nanmean(noStim{ioptoFreq,2}.avgYSpec(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(noStim{ioptoFreq,2}.avgYSpec(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.SpecVC_noStim(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            avgstimPow = nanmean(nanmean(noStim{ioptoFreq,2}.avgPLV(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(noStim{ioptoFreq,2}.avgPLV(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.PLVLPl_VC_noStim(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;

            avgstimPow = nanmean(nanmean(noStim{ioptoFreq,3}.avgPLV(foiMask(iFreq,:),stimtvecMask)));
            avgbasePow = nanmean(nanmean(noStim{ioptoFreq,3}.avgPLV(foiMask(iFreq,:),basetvecMask))); % get 1 value
            sessionAvg.PLVPPC_VC_noStim(ioptoFreq,iFreq,iSession) = (avgstimPow - avgbasePow)/avgbasePow;        
        end
    end
    save([GroupAnalysisDir 'sessionAvg_' freqNames{ioptoFreq} '_20mW.mat'],'sessionAvg', '-v7.3');
end



%%
% plot % power change from baseline


for ioptoFreq = 1:numFreqs
    optoName = freqNames{ioptoFreq};
    colLabel = {'No stim', optoName};
    load([GroupAnalysisDir 'sessionAvg_' optoName '_20mW.mat']);
    
    fig = figure('position', [20    20   600   900]);lw = 2; %x,y,width,height
    yLim = {[-0.3,0.3],[-0.3,1.2],[-0.7,0.8]};
    for iFreq = 1:numFreqs
        subplot(3,3,numFreqs*iFreq-2);
        mat2plot = [squeeze(sessionAvg.SpecLPl_noStim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.SpecLPl(ioptoFreq,iFreq,:))];       
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
        mat2plot = [squeeze(sessionAvg.SpecPPC_noStim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.SpecPPC(ioptoFreq,iFreq,:))];
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
        mat2plot = [squeeze(sessionAvg.SpecVC_noStim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.SpecVC(ioptoFreq,iFreq,:))];
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
        mat2plot = [squeeze(sessionAvg.SpecLPl_noStim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.SpecLPl(ioptoFreq,iFreq,:))];       
        p(numFreqs*iFreq-2) = anova1(mat2plot,colLabel);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
        mat2plot = [squeeze(sessionAvg.SpecPPC_noStim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.SpecPPC(ioptoFreq,iFreq,:))];
        p(numFreqs*iFreq-1) = anova1(mat2plot,colLabel);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
        mat2plot = [squeeze(sessionAvg.SpecVC_noStim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.SpecVC(ioptoFreq,iFreq,:))];
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
        
        mat2plot = [squeeze(sessionAvg.PLVLPl_PPC_noStim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.PLVLPl_PPC(ioptoFreq,iFreq,:))];
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
        mat2plot = [squeeze(sessionAvg.PLVLPl_VC_noStim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.PLVLPl_VC(ioptoFreq,iFreq,:))];
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
        mat2plot = [squeeze(sessionAvg.PLVPPC_VC_noStim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.PLVPPC_VC(ioptoFreq,iFreq,:))];
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
        mat2plot = [squeeze(sessionAvg.PLVLPl_PPC_noStim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.PLVLPl_PPC(ioptoFreq,iFreq,:))];
        p(numFreqs*iFreq-2) = anova1(mat2plot,colLabel);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
        mat2plot = [squeeze(sessionAvg.PLVLPl_VC_noStim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.PLVLPl_VC(ioptoFreq,iFreq,:))];
        p(numFreqs*iFreq-1) = anova1(mat2plot,colLabel);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
        mat2plot = [squeeze(sessionAvg.PLVPPC_VC_noStim(ioptoFreq,iFreq,:)),squeeze(sessionAvg.PLVPPC_VC(ioptoFreq,iFreq,:))];
        p(numFreqs*iFreq) = anova1(mat2plot,colLabel);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
    end
    p_corrected = bonf_holm(p,.05);
    save([GroupAnalysisDir 'sessionAvg_' optoName 'Stim_PLV_20mW_pvalue.mat'],'p','p_corrected');
    clear p p_corrected
    close all
end


% function sub_color_myBars(CT, numDiv)
% a = get(get(gca,'children'),'children');   % Get the handles of all the objects
% t = get(a,'tag');   % List the names of all the objects
% for i=1:numDiv
%     if i == 1 % Set the color of the first box to gray
%         box = a(3*numDiv+1-i); color = [0,0,0]; set(box,'Color',color);
%     else
%         box = a(3*numDiv+1-i);   
%         color = CT(i, :);
%         set(box, 'Color', color);
%     end
% end
% end
