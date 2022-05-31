clear
tic

cluster = 0;
skipRec = 0;
linORlog = 2; %freqs of interest: 1=linear 2=log
MedianorPCA = 3; 
%animalCode = '0168';
animalCode = '0180';

if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
switch animalCode
    case '0168'
        PreprocessDir = ['D:/FerretData/' animalCode '/Preprocessed/'];
        AnalysisDir   = ['D:/FerretData/' animalCode '/Analyzed/'];
        BehavDatDir   = ['D:/FerretData/' animalCode '/behav/'];
        GroupAnalysisDir = ['D:/FerretData/' animalCode '/GroupAnalysis/Spec/'];

    case {'0173','0180','0181'}
        PreprocessDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Preprocessed/'];
        AnalysisDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Analyzed/'];
        BehavDatDir   = ['Z:/Ferret Data/' animalCode '/behav/'];
        GroupAnalysisDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/GroupAnalysis/Spec/'];

end
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/behav/'];
    GroupAnalysisDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/GroupAnalysis/Spec/'];

    %code for initialising parallel computing
    numCore = 24; % USR DEFINE
    myPool = parpool('local',numCore,'SpmdEnabled',false);  

end

% Different folder name if use median or PCA LFP
if MedianorPCA == 0;     folderSuffix = '_validChns'; %use all valid channels
elseif MedianorPCA == 1; folderSuffix = '_median'; %use median of all valid channels
elseif MedianorPCA == 2; folderSuffix = '_PCA';
elseif MedianorPCA == 3; folderSuffix = '_firstChn'; % randomly pick 1 channel
end

fileInfo = dir([PreprocessDir animalCode '_Opto*']); % detect files to load/convert  '_LateralVideo*'

% loop through each recording
for irec = 6%1:numel(fileInfo)
    % so that if one record doesn't work, others can still run
    PulvOpto_rec_cluster(animalCode,irec,fileInfo,folderSuffix, PreprocessDir, AnalysisDir, BehavDatDir, GroupAnalysisDir,...
        cluster, skipRec, linORlog, MedianorPCA);

end % end of all records for an animal

if cluster == 1; delete(myPool); end

x = toc;
fprintf('time required =%f sec\n', x);

