clear all
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

% Directory info
GroupAnalysisDir = 'E:/FerretData/0168/GroupAnalysis/Spec/';

% time info
twin = [-2 2]; %<<<--- interested time window around event
baseTwin = [-1.5 -0.5]; %<<<--- baseline to subtract from spectrogram 
stimTwin = [0 1.2]; %1sec stimulation, 2-3 sec ITI
tvec = is_load('E:\FerretData\0168\GroupAnalysis\Spec\sessions\TrialSpec_001_VC_Left.mat','tvec');
foi  = is_load('E:\FerretData\0168\GroupAnalysis\Spec\sessions\TrialSpec_001_VC_Left.mat','foi');
baseTwins = is_load('E:\FerretData\0168\GroupAnalysis\Spec\sessions\TrialSpec_001_VC_Left.mat','baseTwins');
condIDs = [1 4];
%tvec = twin(1):0.01:twin(2); %make sure tvec length matches avg
basetvecMask = tvec>= baseTwin(1) & tvec<= baseTwin(2);
stimtvecMask = tvec>= stimTwin(1) & tvec<= stimTwin(2);
lowFreq = 2; highFreq = 128; linORlog = 2;
[tickLoc, tickLabel] = getTickLabel(lowFreq, highFreq, 150, linORlog);

% to get gradient color
numDiv     = 5; % divide into 8 equal (in samples) parts
CT=cbrewer('seq', 'Blues', numDiv+2); ColorSet = CT(2:numDiv+1,:); %first 2 colors too light, light->dark
% bm       = blueMap; % total 64 colors
% ColorSet = bm(10:10:10*numDiv,:);
% Opto info
optoFreqs = [4.3, 11.7, 38]; %theta, alpha, gamma
numFreqs  = numel(optoFreqs);
freqNames = {'Theta','Alpha','Gamma'};
optoFtwins = {[3,6],[10,14],[30,60]}; %define freq band of interest
optoAmps  = [10, 20, 30, 40, 50];
numAmps   = numel(optoAmps);
numfoi = length(foi);
numtvec = length(tvec);
numTrials = 50;
for iFreqs = 1:numFreqs
    foiMask(iFreqs,:) = foi>= optoFtwins{iFreqs}(1) & foi<= optoFtwins{iFreqs}(2);
end
%%
%dates with 10-50mW parameterization
sessionDates = {'20180808','20180809','20180814','20180815','20180816','20180817','20180820'};
numDates  = numel(sessionDates);

% Initialize vectors to save
for iRegion = 1:numRegions
    regionName = regionNames{iRegion};
    trialSpecAll.(regionName) = nan(numFreqs, numAmps, numDates*numTrials, numfoi, numtvec);
    trialSpecAllnormed.(regionName) = nan(numFreqs, numAmps, numDates*numTrials, numfoi, numtvec);
end

% loop through each recording
% initialize counter to keep track of how many trials in each condition
trialStartInd = ones(numFreqs, numAmps, numRegions); % to keep track of trial counts
trialEndInd   = zeros(numFreqs, numAmps, numRegions);
trialStartInd_noStim = ones(numFreqs, numAmps, numRegions); 
trialEndInd_noStim   = zeros(numFreqs, numAmps, numRegions);

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
    for ioptoFreq = 1:numFreqs
        for ioptoAmp = 1:numAmps
            sessionID = sessionIDs{ioptoFreq}{ioptoAmp};
            fprintf(['Processing session ' sessionID '\n']);
            for iRegion = 1:numRegions
                regionName = regionNames{iRegion};                
                if strcmp(sessionID, '133'); sessionID = '111';
                elseif strcmp(sessionID, '132'); sessionID = '112';end % fill in 2 missing sessions
                fileName = ['TrialSpec_' sessionID '_' regionName '_Left.mat'];
                spec = is_load([GroupAnalysisDir 'sessions/' fileName],'TrialSpec');
                
                trialEndInd(ioptoFreq,ioptoAmp, iRegion) = trialStartInd(ioptoFreq,ioptoAmp,iRegion) + size(spec.(regionName),1)-1;
                trialSpecAll.(regionName)(ioptoFreq,ioptoAmp,trialStartInd(ioptoFreq,ioptoAmp,iRegion):trialEndInd(ioptoFreq,ioptoAmp,iRegion),:,:) = spec.(regionName);
                trialSpecAllnormed.(regionName)(ioptoFreq,ioptoAmp,trialStartInd(ioptoFreq,ioptoAmp,iRegion):trialEndInd(ioptoFreq,ioptoAmp,iRegion),:,:) = spec.([regionName '_normed']);
                trialStartInd(ioptoFreq,ioptoAmp,iRegion) = trialEndInd(ioptoFreq,ioptoAmp,iRegion) + 1;
                
                if numel(condID) == 2
                fileName = ['TrialSpec_' sessionID '_' regionName '_NoStim.mat'];    
                spec = is_load([GroupAnalysisDir 'sessions/' fileName],'TrialSpec');
                trialEndInd_noStim(ioptoFreq,ioptoAmp,iRegion) = trialStartInd_noStim(ioptoFreq,ioptoAmp,iRegion) + size(spec.(regionName),1)-1;
                trialSpecAll_noStim.(regionName)(ioptoFreq,ioptoAmp,trialStartInd_noStim(ioptoFreq,ioptoAmp,iRegion):trialEndInd_noStim(ioptoFreq,ioptoAmp,iRegion),:,:) = spec.(regionName);
                trialSpecAllnormed_noStim.(regionName)(ioptoFreq,ioptoAmp,trialStartInd_noStim(ioptoFreq,ioptoAmp,iRegion):trialEndInd_noStim(ioptoFreq,ioptoAmp,iRegion),:,:) = spec.([regionName '_normed']);
                trialStartInd_noStim(ioptoFreq,ioptoAmp,iRegion) = trialEndInd_noStim(ioptoFreq,ioptoAmp,iRegion) + 1;
                end
            end
        end
    end
end
% preallocated enough rows, the empty trials were filled with nan, now
% delete those rows in the 3rd dimension using trialEndInd
for iRegion = 1:numRegions
    regionName = regionNames{iRegion};
    fprintf(['Eliminating nan rows and saving for ' regionName '\n']);
    
    trialSpecAll.(regionName)(:,:,trialEndInd(iRegion)+1:end,:,:) = [];
    trialSpecAllnormed.(regionName)(:,:,trialEndInd(iRegion)+1:end,:,:) = [];
    trialSpecAll_noStim.(regionName)(:,:,trialEndInd_noStim(iRegion)+1:end,:,:) = [];
    trialSpecAllnormed_noStim.(regionName)(:,:,trialEndInd_noStim(iRegion)+1:end,:,:) = [];
end
save([GroupAnalysisDir 'trialPool_SpecAll.mat'],'trialSpecAll','trialSpecAll_noStim','-v7.3');
save([GroupAnalysisDir 'trialPool_SpecAllnormed.mat'],'trialSpecAllnormed','trialSpecAllnormed_noStim','-v7.3');
clear spec

%% Calculate and Plot time-resolved baseline normalized power
% for ioptoFreq = 1:numFreqs
%     optoName = freqNames{ioptoFreq};
%     %yLim = {[0.5,1.6],[0.7,2.0],[0.5,2.3]};
%     fig = figure('name',[optoName 'Stim: time-resolved baseline normalized power'],'position', [10 20 370 580]);lw = 2; %x,y,width,height
%     xLim = [-1,2];
%     for iFreq = 1:numFreqs
%         foiName = freqNames{iFreq};
%         for iRegion = 1:numRegions
%             regionName = regionNames{iRegion};
%             subplot(3,3,numFreqs*iFreq-3+iRegion);
%             for iAmp = 1:numAmps
%                 if iAmp == 1 % also plot baseline
%                     noStim2plot = reshape(squeeze(nanmean(trialSpecAllnormed_noStim.(regionName)(:,:,:,foiMask(iFreq,:),:),4)),[],length(tvec));
%                     noStim2plotNormed = pow2db(noStim2plot);
%                     noStim2plotNormed = baselineNorm(noStim2plotNormed, tvec, [-1.5,-0.2],0);
%                     median = nanmedian(noStim2plotNormed,1);
%                     mean = nanmean(noStim2plotNormed,1); % avg over trials make into a row vector
%                     sem  = nanstd(noStim2plotNormed,[],1)/sqrt(size(noStim2plotNormed,1));
%                     noStim2plot_normed(iRegion,iFreq,:,:) = noStim2plotNormed;
%                     shadedErrorBar(tvec, mean, sem, 'k-',0.1); hold on; %last is transparency level            
%                 end
%                     
%                 mat2plot = squeeze(nanmean(trialSpecAllnormed.(regionName)(ioptoFreq,iAmp,:,foiMask(iFreq,:),:),4));
%                 mat2plotNormed = pow2db(mat2plot(1:end-2,:));
%                 mat2plotNormed = baselineNorm(mat2plotNormed,tvec, [-1.5,-0.2],0);
%                 %plot(tvec,mat2plot) % to visualize all trials time course
%                 median   = squeeze(nanmedian(mat2plotNormed,1)); %median across trials
%                 mean     = nanmean(mat2plotNormed,1); %median across trials
%                 sem      = nanstd(mat2plotNormed,[],1)/sqrt(size(mat2plotNormed,1));
%                 Stim2plot_normed(ioptoFreq,iAmp,iRegion,iFreq,1:size(mat2plotNormed,1),:) = mat2plotNormed; % save to a file
%                 Stim2plot_normed(ioptoFreq,iAmp,iRegion,iFreq,size(mat2plotNormed,1):end,:) = nan; % just use pair-wise correlation to avoid nan later
%                 shadedErrorBar(tvec, mean, sem, {'color',ColorSet(iAmp,:)},0.1); hold on; %last is transparency level            
%                 if iFreq == 1; title(regionName); ylim([-2,2]);
%                 elseif iFreq == 2; ylim([-1,2]);
%                 elseif iFreq == 3; ylim([-0.3,1.3]); xlabel('Time [sec]');end
%                 xlim(xLim); 
%                 if iRegion == 1; ylabel(['BLNormed ' foiName ' power']);end                    
%             end % mean looks smoother            
%         end
%         %legend%('10mW','20mW','30mW','40mW','50mW'); 
%     end    
%     savefig(fig, [GroupAnalysisDir 'trialPool_tBLNormedPowerdb_mean_' optoName 'Stim.fig'],'compact');
%     saveas(fig, [GroupAnalysisDir 'trialPool_tBLNormedPowerdb_mean_' optoName 'Stim.png']);
% end
% save([GroupAnalysisDir 'trialPool_tBLNormedPowerdb_mean.mat'],'Stim2plot_normed','noStim2plot_normed','-v7.3');


%% Plot time-resolved power2db and then baseline normalized + sem -- this way baseline is flat
for ioptoFreq = 1:numFreqs
    optoName = freqNames{ioptoFreq};
    %yLim = {[0.5,1.6],[0.7,2.0],[0.5,2.3]};
    fig = figure('name',[optoName 'Stim_trialMean+sem_3Freq_5Amp'], 'position', [10 20 650 900]);lw = 4; %x,y,width,height [10 20 370 580]
    xLim = [-1,2];
    for iFreq = 1:numFreqs
        foiName = freqNames{iFreq};
        for iRegion = 1:numRegions
            regionName = regionNames{iRegion};
            subplot(3,3,numFreqs*iFreq-3+iRegion);
            for iAmp = 1:numAmps
                if iAmp == 1 % also plot baseline
                    noStim2plot = reshape(squeeze(nanmean(trialSpecAll_noStim.(regionName)(:,:,:,foiMask(iFreq,:),:),4)),[],length(tvec));
                    noStim2plotNormed = pow2db(noStim2plot);
                    noStim2plotNormed = baselineNorm(noStim2plotNormed, tvec, [-1.5,-0.2],0);
                    median = nanmedian(noStim2plotNormed,1);
                    mean = nanmean(noStim2plotNormed,1); % avg over trials make into a row vector
                    sem  = nanstd(noStim2plotNormed,[],1)/sqrt(size(noStim2plotNormed,1));
                    noStim2plot_normed(iRegion,iFreq,:,:) = noStim2plotNormed;
                    shadedErrorBar(tvec, mean, sem, {'k-','LineWidth',lw},0.1); hold on; %last is transparency level            
                end
                    
                mat2plot = squeeze(nanmean(trialSpecAll.(regionName)(ioptoFreq,iAmp,:,foiMask(iFreq,:),:),4));
                mat2plotNormed = pow2db(mat2plot(1:end-2,:));
                mat2plotNormed = baselineNorm(mat2plotNormed,tvec, [-1.5,-0.2],0);
                %plot(tvec,mat2plot) % to visualize all trials time course
                median   = squeeze(nanmedian(mat2plotNormed,1)); %median across trials
                mean     = nanmean(mat2plotNormed,1); %median across trials
                sem      = nanstd(mat2plotNormed,[],1)/sqrt(size(mat2plotNormed,1));
                Stim2plot_normed(ioptoFreq,iAmp,iRegion,iFreq,1:size(mat2plotNormed,1),:) = mat2plotNormed; % save to a file
                Stim2plot_normed(ioptoFreq,iAmp,iRegion,iFreq,size(mat2plotNormed,1):end,:) = nan; % just use pair-wise correlation to avoid nan later
                shadedErrorBar(tvec, mean, sem, {'color',ColorSet(iAmp,:),'LineWidth',lw},0.1); hold on; %last is transparency level            
                if iFreq == 1; title(regionName); xticklabels([]); if ioptoFreq ==3 ylim([-2,1]); else ylim([-1,1]); end
                elseif iFreq == 2; ylim([-0.8,1.5]);xticklabels([]);
                elseif iFreq == 3; ylim([-0.3,1.5]); xlabel('Time [sec]');end
                xlim(xLim); 
                if iRegion == 1; ylabel(['BLNormed ' foiName ' power']);else yticklabels([]); end                   
            end % mean looks smoother            
        end
        %legend%('10mW','20mW','30mW','40mW','50mW'); 
    end
    fig.Color = 'white';
    savefig(fig, [GroupAnalysisDir 'trialPool_tBLPowerdb_mean_' optoName 'Stim.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'trialPool_tBLPowerdb_mean_' optoName 'Stim.png']);
end
save([GroupAnalysisDir 'trialPool_tBLPowerdbNormed_mean.mat'],'Stim2plot_normed','noStim2plot_normed','-v7.3');


%% calculate mean theta, alpha, gamma change from baseline using foiMask and stimtvecMask
StimTavg_normed = nanmean(Stim2plot_normed(:,:,:,:,:,stimtvecMask),6); %ioptoFreq x ioptoAmp x iRegion x iFOI x iTrial
noStimTavg_normed = nanmean(noStim2plot_normed(:,:,:,stimtvecMask),4); %iRegion x iFOI x iTrial

save([GroupAnalysisDir 'trialPool_Tavg_normedPowerdb_mean.mat'],'StimTavg_normed','noStimTavg_normed','-v7.3');
% size(StimTavg_normed = 3     5     3     3    300)


%% plot correlation
for ioptoFreq = 1:numFreqs
    optoName = freqNames{ioptoFreq};
    fig = figure('name',[optoName 'Stim: avg power change correlation'],'position', [10 20 370 580]);sz = 2; %x,y,width,height

    for iRegion = 1:numRegions
        regionName = regionNames{iRegion};
        thetaVec = reshape(StimTavg_normed(ioptoFreq,:,iRegion,1,:),numAmps, [])'; % - noStimTavg
        alphaVec = reshape(StimTavg_normed(ioptoFreq,:,iRegion,2,:),numAmps, [])'; % convert to column vector since concatanation is column wise
        gammaVec = reshape(StimTavg_normed(ioptoFreq,:,iRegion,3,:),numAmps, [])';        
        clInd = repmat([1;2;3;4;5], [1,size(thetaVec,1)])'; % same color matrix
        clIndVec = clInd(:)';
        for iTime=1:numel(clIndVec)
            clVec(iTime,:) = ColorSet(clIndVec(iTime),:);
        end
        %keepMask = (~isnan(thetaVec) & ~isnan(alphaVec) & ~isnan(gammaVec);
        % get a linear vector
        subplot(3,3,iRegion)
        scatter(thetaVec(:)', alphaVec(:)',sz, clVec); %thetaVec(:)' gives a row vector with optoAmp 111222333444555
        lsline;
        R = corrcoef(thetaVec(:)',alphaVec(:)','Rows','pairwise'); Rtext = sprintf('%0.2f', R(1,2));
        title([regionName ' corr coef: ' Rtext]);
        xlabel('theta power');
        if iRegion == 1; ylabel(['alpha power']);end

        subplot(3,3,iRegion+numFreqs)
        scatter(thetaVec(:)', gammaVec(:)',sz, clVec); %thetaVec(:)' gives a row vector with optoAmp 111222333444555
        lsline;
        R = corrcoef(thetaVec(:)',gammaVec(:)','Rows','pairwise'); Rtext = sprintf('%0.2f', R(1,2));
        title(['corr coef: ' Rtext]);
        xlabel('theta power');
        if iRegion == 1; ylabel(['gamma power']);end
        
        subplot(3,3,iRegion+2*numFreqs)
        scatter(alphaVec(:)', gammaVec(:)',sz, clVec); %thetaVec(:)' gives a row vector with optoAmp 111222333444555
        lsline;
        R = corrcoef(alphaVec(:)',gammaVec(:)','Rows','pairwise'); Rtext = sprintf('%0.2f', R(1,2));
        title([' corr coef: ' Rtext]);
        xlabel('alpha power');
        if iRegion == 1; ylabel(['gamma power']);end
    end
    savefig(fig, [GroupAnalysisDir 'trialPool_Tavg_Powerdb_change_correlation_' optoName 'Stim.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'trialPool_Tavg_Powerdb_change_correlation_' optoName 'Stim.png']);
end



%% Box plot of T avg
colLabel = {'No stim', '10', '20','30','40','50mW'};
%groupVec = [zeros(1,size(noStimTavg_normed,3)),reshape(repmat([1;2;3;4;5],1,size(StimTavg_normed,5))',1,[])]; %eg. 000111122223333
CTall  = [[0,0,0]; CT(2:6,:)];
%------------- for normalized power
for ioptoFreq = 1:numFreqs
    optoName = freqNames{ioptoFreq};
    fig = figure('name',[optoName 'Stim: trialPool_boxPlot_power_5Amp'], 'position', [10 20 370 580]);lw = 2; %x,y,width,height
    %yLim = {[-2.3,2],[-1,2],[-0.5,1.2]}; %for box plot without whisker
    
    for iRegion = 1:numRegions
        regionName = regionNames{iRegion};
        for iFreq = 1:numFreqs
            subplot(3,3,numFreqs*iFreq-3+iRegion);
            noStimBox = reshape(noStimTavg_normed(iRegion,iFreq,:),1,[]);
            StimBox   = reshape(squeeze(StimTavg_normed(ioptoFreq,:,iRegion,iFreq,:))',1,[]);
            %%StimBox(isnan(StimBox)) = []; %don't need to deal with nan for box plot                   
            mat2boxplot{ioptoFreq,iRegion,iFreq}  = [noStimBox, StimBox];
            boxplot(mat2boxplot{ioptoFreq,iRegion,iFreq},groupVec,'Colors',CTall,'PlotStyle','compact','Whisker',0,'OutlierSize',0.01,'Labels',{'0','10','20','30','40','50'} ); %'BoxStyle','filled',,can't color
            
            % alternativaly, plot mean + sem
            noStimMean = nanmean(noStimTavg_normed(iRegion,iFreq,:),3);
            noStimSem  = nanstd(noStimTavg_normed(iRegion,iFreq,:),[],3)/sqrt(size(noStimTavg_normed,3));
            StimMean   = squeeze(nanmean(StimTavg_normed(ioptoFreq,:,iRegion,iFreq,:),5)');
            StimSem    = squeeze(nanstd(StimTavg_normed(ioptoFreq,:,iRegion,iFreq,:),[],5)');
            errorbar([noStimMean StimMean],[noStimSem StimSem],CTAll);
            
            %color_myBars(CT, numDiv+1);
            ylim(yLim{iFreq}); 
            ylabel([freqNames{iFreq} ' power']); % add ylabel distorts figure shape
            if iFreq == 1; title(regionName); %legend('no stim','10mW','20mW','30mW','40mW','50mW'); 
            end
        if iFreq == 3; xlabel('Amplitude [mW]');end % add xlabel distorts figure shape
        end
    end
    savefig(fig, [GroupAnalysisDir 'trialPool_boxPlot_' optoName 'Stim_power_5Amp_zoom.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'trialPool_boxPlot_' optoName 'Stim_power_5Amp_zoom.png']); 
end
save([GroupAnalysisDir 'trialPool_boxPlot_power_5Amp.mat'],'mat2boxplot', '-v7.3');


for ioptoFreq = 1:numFreqs
    % one-way anova comparing no stim and stim, save p values for correcting multiple comparison
    %p = 3 optoFreq x (3 region x 3 freqs) x (1 main p + 5 pairwise p)   = 3x9x6  
    for iRegion = 1:numRegions
        for iFreq = 1:numFreqs
            [p(ioptoFreq,iRegion,iFreq),~,stats] = anova1(mat2boxplot{ioptoFreq,iRegion,iFreq},groupVec);%produce anova table and a boxplot with notch on (same as boxplot(y,labels,'notch','on')
            [c,~,~,gnames] = multcompare(stats);
            p(ioptoFreq,iRegion,iFreq,2:6) = c(1:5,6)'; %7th column is p value, 1-5 row is no stim vs. 10-50mA
            p_corrected(ioptoFreq,iRegion,iFreq,:) = bonf_holm(squeeze(p(ioptoFreq,iRegion,iFreq,2:6)),.05); % only correct for pair-wise comparison
        end
    end
end
save([GroupAnalysisDir 'trialPool_boxPlot_power_5Amp_pvalue.mat'],'p','p_corrected');
clear p p_corrected



function sub_color_myBars(CT, numDiv)
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects
for i=1:numDiv
    if i == 1 % Set the color of the first box to gray
        box = a(3*numDiv+1-i); color = [0,0,0]; set(box,'Color',color);
    else
        box = a(3*numDiv+1-i);   
        color = CT(i, :);
        set(box, 'Color', color);
    end
end
end
