function s1_preprocessEEG(subject)
% s1_preprocessEEG
% Plugins required: ERPLAB, ICLabel, Viewprops, icablinkmetrics
%
% USAGE EXAMPLES
% preprocessEEG('001')
%
% The function takes only 1 input argument (character; e.g., '001').
%
% Written in MATLAB R2018b
% Last modified by Hause Lin 09-05-19 8:51 PM hauselin@gmail.com

%% Specify paths

% Change the paths below accordingly before you use this script!
PATHS = struct();
PATHS.cwd = pwd;
PATHS.rawdatadir = '/Users/hause/Dropbox/Working Projects/Akina Physical Effort/Zaibas RewP study/DataRawMat/';
PATHS.outputdir = '/Users/hause/Dropbox/Working Projects/Akina Physical Effort/Zaibas RewP study/DataPreprocessed/';
PATHS.log = fullfile(PATHS.cwd,'Log'); % creates a Log directory in your current directory

% addpath to eeglab (in case)
% addpath('/psyhome/u4/linhause/matlabtoolboxes/eeglab14_1_2b'); % addpath on UTSC cluster
% addpath('/users/hause/dropbox/Apps/MATLAB Toolboxes and Packages/eeglab14_1_2b') % add path on local machine

%% Specify preprocessing parameters

% Change the parameters for your preprocessing pipeline!
cfg = struct();

% filter settings
cfg.highPassForICA = 1; % high pass filter for ICA training only
cfg.highPass = 0.1; % actual high pass (for use after ICA decomposition)
% cfg.lowPass = 40; % usually for ICA training only (maybe unnecessary anyway)

% trim/remove long sections of data without events (e.g., between-block breaks)
cfg.minTimeThreshold = 5000; % minimum time (s) between two events
cfg.bufferTime = 4500; % add buffer time (s) (pre/post event)

cfg.epochArtifactThresholdZValue = 4; % SD cut-off for pre-ICA epoch artifact rejection 
cfg.badChanThresholdSD = 5; % SD cutoff for bad electrode detection

cfg.eyeChan = {'veog1','veog2','heogL','heogR'}; % name of eye/ocular channels in data
cfg.eyeChanBESA = {'SO2','IO2','LO1','LO2'}; % name of eye channels (BESA names) (make sure order matches cfg.eyeChan above!)
cfg.emgChan = {'CorsOut','CorsIns','CORRins','CORRout','ZYGup','ZYGlow','COORins','COORout'}; % name of emg channels (to be removed in this preprocessing pipeline)

% dummy epoch duration for ICA (1 second is good!)
cfg.epochDuration = 1; % epochs data into 1s windows to detect artifacts (and remove epochs) before ICA

%% Perform checks before preprocessing

% If everything is set up properly, you shouldn't have to edit any script
% from here onwards.

if ~exist(PATHS.outputdir) 
    mkdir(PATHS.outputdir);
end

if ~exist(PATHS.log) 
    mkdir(PATHS.log);
end

filetoread = dir(fullfile(PATHS.rawdatadir,['*_' subject '.mat']));
if length(filetoread) ~= 1 % if more or less than 1 matching file found, skip
    dlmwrite(fullfile(PATHS.log,['Error_check_rawData_' subject '.txt']), ['Subject ' subject '. Incorrect number of files in raw data directory.'], 'delimiter','');
    return
end

disp(['Preprocessing subject ' subject, ' with ICA now...']);

%% Start preprocessing

try

%% load data    

load(fullfile(PATHS.rawdatadir,filetoread.name)); % read data    
EEG.subject = subject;

%% Edit channels

allChannels = {EEG.chanlocs.labels}; % all channel labels

if ~isempty(cfg.eyeChan) % edit/rename eye/ocular channels to match BESA channel names
    for chanI = 1:length(cfg.eyeChan)
        eyeChanIdx = find(strcmpi(cfg.eyeChan{chanI}, allChannels)); % get index of eye channel
        if isempty(eyeChanIdx) % if no index found, skip to next iteration
            % do nothing
        else
            EEG = pop_chanedit(EEG,'changefield',{eyeChanIdx 'labels' cfg.eyeChanBESA{chanI}}); % rename channels to BESA channel names
        end
    end
end
    
% remove emg channels
if ~isempty(cfg.emgChan)
    toRemove = find(ismember(allChannels,cfg.emgChan));
    toRemove = allChannels(toRemove); % channel labels to remove
    EEG = pop_select(EEG,'nochannel',toRemove); % remove EMG channels
    EEG = eeg_checkset(EEG);
end
EEG.comments = pop_comments(EEG.comments,'','Removed EMG channels.',1);

% add channel locations
EEG = pop_chanedit(EEG,'lookup','standard-10-5-cap385.elp');
EEG.comments = pop_comments(EEG.comments,'','Added channel locations.',1);

% make a copy of EEG data structure for post ICA
EEGcopy = EEG;

%% Band pass filter data

% highpass filter for ICA (1 Hz)
EEG = pop_basicfilter(EEG,1:length(EEG.chanlocs),'Boundary','boundary','Cutoff',cfg.highPassForICA,'Design','butter','Filter','highpass','Order',2,'RemoveDC','on');

% or bandpass filter for ICA?
% EEG = pop_basicfilter(EEG,1:length(EEG.chanlocs),'Boundary','boundary','Cutoff',[cfg.highPassForICA cfg.lowPass],'Design','butter','Filter','bandpass','Order',2,'RemoveDC','on');

% highpass filter for actual analyses (0.1 Hz)
EEGcopy = pop_basicfilter(EEGcopy,1:length(EEGcopy.chanlocs),'Boundary','boundary','Cutoff',cfg.highPass,'Design','butter','Filter','highpass','Order',2,'RemoveDC','on');

%% Trim data using ERPLAB function

EEG = pop_erplabDeleteTimeSegments(EEG,'displayEEG',0,'endEventcodeBufferMS',cfg.bufferTime,'ignoreUseType','ignore','startEventcodeBufferMS',cfg.bufferTime,'timeThresholdMS',cfg.minTimeThreshold);
EEGcopy = pop_erplabDeleteTimeSegments(EEGcopy,'displayEEG',0,'endEventcodeBufferMS',cfg.bufferTime,'ignoreUseType','ignore','startEventcodeBufferMS',cfg.bufferTime,'timeThresholdMS',cfg.minTimeThreshold);

%% Convert triggers/events from string to numeric

% pop_squeezevents(EEG); % summarize events (ERPLAB)
% eeg_eventtypes(EEG) % summarize events (EEGLAB)

EEG = pop_creabasiceventlist(EEG,'AlphanumericCleaning','on','BoundaryNumeric',{-99},'BoundaryString',{'boundary'});
EEGcopy = pop_creabasiceventlist(EEGcopy,'AlphanumericCleaning','on','BoundaryNumeric',{-99},'BoundaryString',{'boundary'});

%% Create dummy epochs for ICA 

% add dummy trigger for epoching
for timeI = 1:(EEG.srate * cfg.epochDuration):EEG.pnts % for every epochDuration (in seconds), add dummy event
     EEG.event(end+1).type = 9999; % add event type
     EEG.event(end).latency = timeI; % add event latency for the new event added above
end
EEG = eeg_checkset(EEG);
pop_squeezevents(EEG); % summarize events (ERPLAB)

% epoch data
EEG = pop_epoch(EEG,{9999},[0 cfg.epochDuration],'newname','temporary epochs','epochinfo', 'yes');

%% Reject epochs with artifact before ICA

% identify which channels to look for artifacts in
allChannels = {EEG.chanlocs.labels}; % all channel labels
channelsForArtifactDetection = find(~ismember(allChannels, {'SO1' 'SO2' 'IO1' 'IO2' 'LO1' 'LO2' 'Fp1' 'Fpz' 'Fp2'})); % exclude these channels

% identify which epochs contain artifact
try
    EEG = pop_jointprob(EEG,1,channelsForArtifactDetection,cfg.epochArtifactThresholdZValue,cfg.epochArtifactThresholdZValue,0,0);
end
try
    EEG = pop_rejkurt(EEG,1,channelsForArtifactDetection,cfg.epochArtifactThresholdZValue,cfg.epochArtifactThresholdZValue,0,0);
end
EEG = eeg_rejsuperpose(EEG,1,1,1,1,1,1,1,1); % update EEG.reject.rejglobal and EEG.reject.rejglobalE fields with all rejected epochs

% find(EEG.reject.rejglobal) % epochs rejected
% length(find(EEG.reject.rejglobal)) % epochs rejected
% percent epochs rejected
% epochsDeletedPercent = sum(EEG.reject.rejglobal) / length(EEG.reject.rejglobal) * 100;

%% Reject bad epochs before ICA

EEG = pop_rejepoch(EEG,EEG.reject.rejglobal,0);

%% Reject bad channels before ICA

try 
    allChannels = {EEG.chanlocs.labels}; % all channel labels
    channelsForArtifactDetection = find(~ismember(allChannels, {'SO1' 'SO2' 'IO1' 'IO2' 'LO1' 'LO2'})); % exclude these channels when detecting bad electrode
    indelec = 0; % index of bad electrode
    [EEG,indelec,measure] = pop_rejchan(EEG,'threshold',cfg.badChanThresholdSD,'measure','prob','norm','on','elec',channelsForArtifactDetection); % reject bad channels based on probability
    if indelec ~= 0 % if bad channel detected, save information
        badChannelIdx = channelsForArtifactDetection(indelec); % bad channel index
        badChannelName = allChannels(badChannelIdx); % bad channel name
        EEG.badchannel = badChannelName;
        EEG.comments = pop_comments(EEG.comments,'','Removed bad channels.',1);
    end
catch
    indelec = 0; % index of bad electrode
end

%% Run ICA

% EEG = pop_runica(EEG,'extended',1,'icatype','runica'); % runica algorithm
EEG = pop_runica(EEG,'icatype','jader','dataset',1);  % jade algorithm
EEG = eeg_checkset(EEG);

% save/transfer ICA weights
EEGcopy.icawinv = EEG.icawinv;
EEGcopy.icasphere = EEG.icasphere;
EEGcopy.icaweights = EEG.icaweights;
EEGcopy.icachansind = EEG.icachansind;
% EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:); % recompute ICA component channel activation/time-series

EEG = EEGcopy;
clear EEGcopy

%% Remove bad channels in ICA-ed data

% remove bad channel if there's any
if indelec ~= 0
    EEG = pop_select(EEG,'nochannel',badChannelName);
    EEG.comments = pop_comments(EEG.comments,'','Removed bad channels.',1);
    EEG.badchannel = badChannelName;
end
EEG.comments = pop_comments(EEG.comments,'','Performed ICA.',1);

%% Label and mark artifact IC components

try % use iclabel to label components
    EEG = iclabel(EEG); % label ICA components
end
try % use icablinkmetrics to identify blink components
    EEG.icaquant = icablinkmetrics(EEG,'ArtifactChannel',EEG.data(find(strcmp({EEG.chanlocs.labels},'SO2')),:),'Alpha',0.001, 'VisualizeData','False');
end
EEG.comments = pop_comments(EEG.comments,'','Marked artifact IC components',1);
EEG = eeg_checkset(EEG);

%% Save dataset

EEG.history = [];
EEG.setname = [subject '_continuous_icaed'];
EEG.cfg_preprocess = cfg;
save(fullfile(PATHS.outputdir,[EEG.setname '.mat']),'EEG');
dlmwrite(fullfile(PATHS.log,['Finished_Subject_' subject '.txt']),'Finished!','delimiter','');

catch ME % if errors when running script, catch them

save(fullfile(PATHS.log,['Error_MException_' subject '.mat']),'ME'); % save MException object
    
end % end try catch

clear; close all; clc;

end % end  preprocessEEG function