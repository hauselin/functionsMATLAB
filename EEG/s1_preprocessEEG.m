function s1_preprocessEEG(subject)
% Preprocesses EEG data.
% Plugins used: ERPLAB, ANTeepimport (to read eeprobe .cnt)
% Plugins for ICA: ADJUST, ICLabel, Viewprops, icablinkmetrics
%
% USAGE EXAMPLES
% preprocessEEG('001')
%
% The function takes only 1 input argument (character; e.g., '001').
% Assumes that you have a subdirectory called EEGData (can be
% changed below). This EEGData subdirectory should contain one directory per
% subject, and each subject's directory name should match exactly
% the input argument provided to this function. This function
% tries to read raw data from the 'raw' directory.
%
% These folders are required within each subject's directory
% (e.g., subject 001)
% - directory with raw data: ./EEGData/001/raw/
% - directory for output data: ./EEGData/001/continuous/
% - directory for processing parameters: ./EEGData/001/parameters/
%
% MORE EXAMPLES
% 1. Example directory structure: './EEGData/001/raw/Subject001.cnt'
% RUN: preprocessEEG('001')
% 2. Example directory structure: './EEGData/012/raw/Subject012.cnt'
% RUN: preprocessEEG('012')
% 3. Example directory structure: './EEGData/020/raw/Subject20.cnt'
% Note the difference: '020' and 'Subject20.cnt'. This function
% ignores the actual file name ('Subject20.cnt') but instead
% tries to match '020'.
% RUN: preprocessEEG('020') % not preprocessEEG('20')
%
% If useICA = true (default true), ICA decomposition will be run on 1-40Hz
% bandpass filtered data. The ICA weights will then be applied to
% 0.1Hz highpass filtered data. If you don't want to use ICA,
% make sure to set useICA = false.
%
% Preprocessing steps
% More to come...
%
% Written in MATLAB R2016b (also works with 2017b)
% Last modified by Hause Lin 18-05-18 23:02 hauselin@gmail.com

%% Specify variables to set up workspace and preprocessing parameters

% run ICA?
useICA = false;

% keep existing already preprocessed output files in continuous directory?
keepPreprocessed = false; % (default true: keep/don't overwite existing files; false: run preprocessing again and overwrite existing files

% stimulus ttls/event codes
stimulusEvents = [11 12]; % first stimulus of each trial if each trial has multiple stimuli
otherEventCodes = []; % optional; document other event codes (optional) but good to remember what's what

eegdataformat = '.cnt'; % raw data file format
dataDirectory = 'EEGData'; % directory containing subject/participant data

% filter settings
lowPass = 40; % for ICA training only
highPass = 0.1; % actual high pass (for use after ICA decomposition)
highPassForICA = 1; % high pass filter for ICA training only

% trim/remove long sections of data without events (e.g., between-block breaks)
minTimeThreshold = 5000; % minimum time (s) between two events
bufferTime = 4500; % add buffer time (s) (pre/post event)

% SD cut-off for pre-ICA artifact rejection
artifactThresholdZValue = 4;

% dummy epoching for ICA (epoch duration in seconds)
epochDuration = 1; % epochs data into 1s windows to detect artifacts (and remove epochs) before ICA

% emg channel (to be removed from this preprocessing)
emgChan = {'CorsOut','CorsIns','CORRins','CORRout','ZYGup','ZYGlow','COORins','COORout'}; % name of emg channels

% eye/ocular channels
eyeChan = {'veog1','veog2','heogL','heogR'}; % name of eye channels
eyeChanBESA = {'SO2','IO2','LO1','LO2'}; % name of eye channels (BESA names)

% addpath to eeglab (in case)
addpath('/psyhome/u4/linhause/matlabtoolboxes/eeglab14_1_2b'); % addpath on UTSC cluster
addpath('/users/hause/dropbox/Apps/MATLAB Toolboxes and Packages/eeglab14_1_2b') % add path on local machine

%% Check if subject needs to be pre-processed

% Do not modify from here onwards unless you know what you're doing
dataDirectoryAbsPath = fullfile(pwd, dataDirectory);

clc;
if ~exist(fullfile(dataDirectoryAbsPath, subject),'dir') % if directory folder doesn't exist, skip this subject
    disp(['Subject ' subject, ' directory is unavailable']);
    return
end
disp(['Subject ' subject, '']);

if keepPreprocessed % if keeping existing files (don't overwrite)

    filesInContinuousDirectory = dir(fullfile(dataDirectoryAbsPath, subject, 'continuous'));
    match1_ica = sum(~cellfun(@isempty, strfind({filesInContinuousDirectory.name}, '_continuous_ICAweights.set'))) > 0;
    match2_noica = sum(~cellfun(@isempty, strfind({filesInContinuousDirectory.name}, '_continuous.set'))) > 0;

    if useICA && match1_ica % if use ICA
        disp(['Subject ' subject, ' already preprocessed with ICA. Skipping this subject...']);
        dlmwrite(['Skipped_Subject_' subject '_ICA.txt'], 'Skipped','delimiter','');
        return % escape function
    elseif ~useICA && match2_noica
        disp(['Subject ' subject, ' already preprocessed without ICA. Skipping this subject...']);
        dlmwrite(['Skipped_Subject_' subject '_noICA.txt'], 'Skipped','delimiter','');
        return % escape function
    else
        if useICA
           disp(['Preprocessing subject ' subject, ' with ICA now...']);
        else
           disp(['Preprocessing subject ' subject, ' without ICA now...']);
        end
    end

else
    if useICA
       disp(['Preprocessing subject ' subject, ' with ICA now...']);
    else
       disp(['Preprocessing subject ' subject, ' without ICA now...']);
    end
end

%% Start preprocessing

try

    % find files to load/read
    currentSubjectDirectory = fullfile(dataDirectoryAbsPath, subject);
    rawDataDirectory = fullfile(currentSubjectDirectory,'raw'); % directory with raw data
    filesInRawDirectory = dir(rawDataDirectory); % files in raw data directory
    fileIdx = find(~cellfun(@isempty, strfind({filesInRawDirectory.name}, eegdataformat))); % index of file to read

    if length(fileIdx) ~= 1 % if more than one raw data found, skip to next person!
        % save error message to directory
        disp('incorrect number of files in raw data directory');
        tempTxt = ['Subject ' subject '. Incorrect number of files in raw data directory.'];
        dlmwrite(['Error_check_rawData_' subject '.txt'], tempTxt, 'delimiter','');
        % continue % passes control to next iteraction of loop
        return % quit function
    end

    rawFileAbsDirectory = fullfile(rawDataDirectory, filesInRawDirectory(fileIdx).name); % absolute directory path to raw data

    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % run eeglab
    % EEG = pop_loadeep(rawFileAbsDirectory,'triggerfile','on'); % read data into MATLAB/eeglab
    EEG = pop_loadeep_v4(rawFileAbsDirectory); % read data into MATLAB/eeglab
    EEG.subject = subject;

    % pop_squeezevents(EEG); % summarize events (ERPLAB)
    % eeg_eventtypes(EEG) % summarize events (EEGLAB)

    %% Edit channels

    allChannels = {EEG.chanlocs.labels}; % all channel labels

    % edit/rename eye/ocular channels to match BESA channel names
    if ~isempty(eyeChan)
        for chanI = 1:length(eyeChan)
            eyeChanIdx = find(strcmpi(eyeChan{chanI}, allChannels)); % get index of eye channel
            if isempty(eyeChanIdx) % if no index found, skip to next iteration
                % do nothing
            else
                EEG = pop_chanedit(EEG,'changefield',{eyeChanIdx 'labels' eyeChanBESA{chanI}}); % rename channels to BESA channel names
            end
        end
    end
    
    % remove emg channels
    if ~isempty(emgChan)
        toRemove = [];
        for emgChanI = emgChan
            toRemove = [toRemove find(strcmpi(emgChanI{1}, allChannels))]; % find index of emg channels
        end
        toRemove = allChannels(toRemove); % channel labels to remove
        if ~isempty(toRemove)
            EEG_emg = pop_select(EEG,'channel',toRemove); % save EMG channels
        end
        EEG = pop_select(EEG,'nochannel',toRemove); % remove EMG channels
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG,EEG,1,'overwrite','on'); % modify ALLEEG structure, check consistency
    end
    EEG.comments = pop_comments(EEG.comments,'','Removed EMG channels.',1);

    % add channel locations
    EEG = pop_chanedit(EEG,'lookup','standard-10-5-cap385.elp');
    EEG.comments = pop_comments(EEG.comments,'','Added channel locations.',1);

    %% if doing ICA, run the following

    if useICA

        %% Band pass filter for ICA (ERPLAB filter)

        EEG = pop_basicfilter(EEG,1:length(EEG.chanlocs),'Boundary','boundary','Cutoff',[highPassForICA lowPass],'Design','butter','Filter','bandpass','Order', 2,'RemoveDC', 'on');
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,1); % modify ALLEEG structure; check consistency; store in dataset1

        %% Trim data (ERPLAB function)

        EEG = pop_erplabDeleteTimeSegments(EEG,'displayEEG',0,'endEventcodeBufferMS',bufferTime,'ignoreUseType','ignore','startEventcodeBufferMS',bufferTime,'timeThresholdMS',minTimeThreshold);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,1); % modify ALLEEG structure; check consistency; store in dataset1
        
        %% Checking events

        % pop_squeezevents(EEG); % summarize events (ERPLAB)
        % eeg_eventtypes(EEG) % summarize events (EEGLAB)
        % convert triggers from string to numeric
        EEG = pop_creabasiceventlist(EEG,'AlphanumericCleaning','on','BoundaryNumeric',{-99},'BoundaryString',{'boundary'});

        %% Create dummy epochs

        % add dummy trigger for epoching
        for timeI = 1:(EEG.srate * epochDuration):EEG.pnts % for every epochDuration (in seconds), add dummy event
             EEG.event(end+1).type = 9999; % add event type
             EEG.event(end).latency = timeI; % add event latency for the new event added above
        end
        EEG = eeg_checkset(EEG);
        pop_squeezevents(EEG); % summarize events (ERPLAB)

        % epoch data
        EEG = pop_epoch(EEG,{9999},[0 epochDuration],'newname','temporary epochs','epochinfo', 'yes');

        %% Reject epochs with artifact

        % identify which channels to look for artifacts in
        allChannels = {EEG.chanlocs.labels}; % all channel labels
        channelsForArtifactDetection = find(~ismember(allChannels, {'SO1' 'SO2' 'IO1' 'IO2' 'LO1' 'LO2' 'Fp1' 'Fpz' 'Fp2'})); % exclude these channels

        % identify which epochs contain artifact
        try
            EEG = pop_jointprob(EEG,1,channelsForArtifactDetection,artifactThresholdZValue,artifactThresholdZValue,0,0);
        end
        try
            EEG = pop_rejkurt(EEG,1,channelsForArtifactDetection,artifactThresholdZValue,artifactThresholdZValue,0,0);
        end
        EEG = eeg_rejsuperpose(EEG,1,1,1,1,1,1,1,1); % update EEG.reject.rejglobal and EEG.reject.rejglobalE fields with all rejected epochs

        % find(EEG.reject.rejglobal) % epochs rejected
        % length(find(EEG.reject.rejglobal)) % epochs rejected

        % percent epochs rejected
        epochsDeletedPercent = sum(EEG.reject.rejglobal) / length(EEG.reject.rejglobal) * 100;
        subjectVariable = {subject};
        % save as csv in the parameters folder
        writetable(table(subjectVariable , epochsDeletedPercent), fullfile(currentSubjectDirectory, 'parameters', [subject '_preICA_Info.csv']))

        %% Save trial and event information (for debugging purposes)

        epochsRejected = find(EEG.reject.rejglobal); % epochs rejected
        eventsRejectedIdx = find(ismember([EEG.event.epoch], epochsRejected)); % indices of rejected events (each rejected epoch can contain many events)
        eventsRejectedStruct = EEG.event(eventsRejectedIdx); % events deleted structure
        eventsRejectedStructWithout9999 = eventsRejectedStruct([eventsRejectedStruct.type] ~= 9999); % exclude 9999

        eventTable = struct2table(EEG.EVENTLIST.eventinfo); % convert EEG/ERPLAB raw event list to table
        rowIdxToDisable = find(ismember([EEG.EVENTLIST.eventinfo.spoint], [eventsRejectedStructWithout9999.latency])); % which rows to set enable = 0
        eventTable{rowIdxToDisable, 'enable'} = 0;

        % add trial number to table
        eventTable{:, 'trialAll'} = NaN;
        stimulusOnsetIdx = find(ismember(eventTable.code, stimulusEvents));
        nTrials = length(stimulusOnsetIdx);
        eventTable{stimulusOnsetIdx, 'trialAll'} = [1:nTrials]';

        % fill NaN with value from previous index
        replaced = eventTable.trialAll; % make a copy
        for missingI = 1:length(replaced) % for each value in b
            if missingI == 1 % skip first number
                continue
            elseif isnan(replaced(missingI)) % if is NaN, replace with value from previous index
                previousValue = replaced (missingI-1);
                replaced(missingI) = previousValue;
            else
                continue
            end
        end

        eventTable.trialAll = replaced; % replace variable with new variable
        clear replaced;

        eventTable{:, 'trialRemain'} = eventTable.trialAll;
        eventTable{eventTable.enable == 0, 'trialRemain'} = NaN;
        eventTable.codelabel = []; eventTable.binlabel = []; eventTable.dura = []; eventTable.flag = []; eventTable.bini = []; eventTable.bepoch = [];

        % save as csv in the parameters folder
        writetable(eventTable, fullfile(currentSubjectDirectory, 'parameters', [subject '_preICA_eventlist.csv']))

        %% Reject bad epochs

        EEG = pop_rejepoch(EEG,EEG.reject.rejglobal,0);

    end

    %% Reject bad channels

    try
        % find bad channel(s) (5 SD)
        allChannels = {EEG.chanlocs.labels}; % all channel labels
        channelsForArtifactDetection = find(~ismember(allChannels, {'SO1' 'SO2' 'IO1' 'IO2' 'LO1' 'LO2'})); % exclude these channels when detecting bad electrode

        indelec = 0; % index of bad electrode
        thresholdSD = 5; % SD cutoff for bad electrode detection
        [EEG,indelec,measure] = pop_rejchan(EEG,'threshold',thresholdSD,'measure','prob','norm','on','elec',channelsForArtifactDetection); % reject bad channels based on probability
        if indelec ~= 0 % if bad channel detected, save information
            badChannelIdx = channelsForArtifactDetection(indelec); % bad channel index
            badChannelName = allChannels(badChannelIdx); % bad channel name
            writetable(table(subjectVariable, epochsDeletedPercent, badChannelIdx), fullfile(currentSubjectDirectory, 'parameters', [subject '_preICA_Info.csv'])) % save info as csv
            EEG.badchannel = badChannelName;
            EEG.comments = pop_comments(EEG.comments,'','Removed bad channels.',1);
        end
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG,EEG,CURRENTSET,'overwrite','on'); % modify ALLEEG structure, check consistency
    catch
        indelec = 0; % index of bad electrode
    end

    if useICA

        %% Run ICA to get weights

        EEG = pop_runica(EEG,'extended',1,'icatype','runica');
        EEG = eeg_checkset(EEG);

        % save ICA weights
        tempICAweights.icawinv = EEG.icawinv;
        tempICAweights.icasphere = EEG.icasphere;
        tempICAweights.icaweights = EEG.icaweights;
        tempICAweights.icachansind = EEG.icachansind;
        save(fullfile(currentSubjectDirectory, 'parameters', [subject '_ICAweights.mat']), 'tempICAweights') % save ICA weights in parameters directory

        %% Reload data and clean it again with new highpass filter later on

        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % run eeglab
        % EEG = pop_loadeep(rawFileAbsDirectory,'triggerfile','on'); % read data into MATLAB/eeglab
        EEG = pop_loadeep_v4(rawFileAbsDirectory); % read data into MATLAB/eeglab
        EEG.subject = subject;

        %% Edit channels

        allChannels = {EEG.chanlocs.labels}; % all channel labels

        % edit/rename eye/ocular channels to match BESA channel names
        if ~isempty(eyeChan)
            for chanI = 1:length(eyeChan)
                eyeChanIdx = find(strcmpi(eyeChan{chanI}, allChannels)); % get index of eye channel
                if isempty(eyeChanIdx) % if no index found, skip to next iteration
                    % do nothing
                else
                    EEG = pop_chanedit(EEG,'changefield',{eyeChanIdx 'labels' eyeChanBESA{chanI}}); % rename channels to BESA channel names
                end
            end
        end

        % remove emg channels
        if ~isempty(emgChan)
            toRemove = [];
            for emgChanI = emgChan
                toRemove = [toRemove find(strcmpi(emgChanI{1}, allChannels))]; % find index of emg channels
            end
            toRemove = allChannels(toRemove); % channel labels to remove
            if ~isempty(toRemove)
                EEG_emg = pop_select(EEG,'channel',toRemove); % save EMG channels
            end
            EEG = pop_select(EEG,'nochannel',toRemove); % remove EMG channels
            EEG.comments = pop_comments(EEG.comments,'','Removed EMG channels.',1);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG,EEG,CURRENTSET,'overwrite','on'); % modify ALLEEG structure, check consistency
        end

        % remove bad channel if there's any
        if indelec ~= 0
            EEG = pop_select(EEG,'nochannel',badChannelName);
            EEG.comments = pop_comments(EEG.comments,'','Removed bad channels.',1);
            EEG.badchannel = badChannelName;
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG,EEG,CURRENTSET,'overwrite','on'); % modify ALLEEG structure, check consistency
        end

        % add channel location
        EEG = pop_chanedit(EEG,'lookup','standard-10-5-cap385.elp');
        EEG.comments = pop_comments(EEG.comments,'','Added channel locations.',1);

    end

    %% Band pass filter (only high pass at 0.1; no low pass)

    EEG = pop_basicfilter(EEG,1:length(EEG.chanlocs),'Boundary','boundary','Cutoff',highPass,'Design','butter','Filter','highpass','Order',2,'RemoveDC','on');
    EEG.comments = pop_comments(EEG.comments,'',['Highpass filtered at ' num2str(highPass) 'Hz.'],1);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,1); % modify ALLEEG structure; check consistency; store in dataset1

    %% Trim data (ERPLAB function)

    % trim EEG data
    EEG = pop_erplabDeleteTimeSegments(EEG,'displayEEG',0,'endEventcodeBufferMS',bufferTime,'ignoreUseType','ignore','startEventcodeBufferMS',bufferTime,'timeThresholdMS',minTimeThreshold);
    EEG.comments = pop_comments(EEG.comments,'','Trimmed data by deleting eventless time segments.',1);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,1); % modify ALLEEG structure; check consistency; store in dataset1
    
    % trim EMG data
    if ~isempty(toRemove)
        EEG_emg = pop_erplabDeleteTimeSegments(EEG_emg,'displayEEG',0,'endEventcodeBufferMS',bufferTime,'ignoreUseType','ignore','startEventcodeBufferMS',bufferTime,'timeThresholdMS',minTimeThreshold);
    end
    
    %% Checking events

    % pop_squeezevents(EEG); % summarize events (ERPLAB)
    % eeg_eventtypes(EEG) % summarize events (EEGLAB)
    % convert triggers from string to numeric
    EEG = pop_creabasiceventlist(EEG,'AlphanumericCleaning','on','BoundaryNumeric',{-99},'BoundaryString',{'boundary'});
    if ~isempty(toRemove)
        EEG_emg = pop_creabasiceventlist(EEG_emg,'AlphanumericCleaning','on','BoundaryNumeric',{-99},'BoundaryString',{'boundary'});
    end
    
    %% Apply ICA weights from 1-40Hz filtered data to 0.1 Hz filtered data

    if useICA

        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG.icawinv = tempICAweights.icawinv;
        EEG.icasphere = tempICAweights.icasphere;
        EEG.icaweights = tempICAweights.icaweights;
        EEG.icachansind = tempICAweights.icachansind;

        EEG.comments = pop_comments(EEG.comments,'','Finished ICA.',1);
        EEG = eeg_checkset(EEG);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

        % EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:); % recompute ICA component channel activation/time-series

    %% Label and mark artifact IC components

        % if more than 20 electrodes, use ADJUST to find artifact components
        if length(EEG.chanlocs) > 20
            try
                EEG = interface_ADJ(EEG,fullfile(currentSubjectDirectory, 'parameters', [EEG.subject '_adjustICA_report.txt']),0); % don't plot
            catch
                EEG = interface_ADJ(EEG,fullfile(currentSubjectDirectory, 'parameters', [EEG.subject '_adjustICA_report.txt']));
            end
        end

        % use iclabel to label components
        try
            EEG = iclabel(EEG); % label ICA components
        catch
            disp('Failed to label ICA components with iclabel');
        end

        % use icablinkmetrics to identify blink components
        try
            EEG.icaquant = icablinkmetrics(EEG,'ArtifactChannel',EEG.data(find(strcmp({EEG.chanlocs.labels},'SO2')),:),'Alpha',0.001, 'VisualizeData','False');
        catch
            disp('Failed to identify blink components with icablinkmetrics');
        end

        EEG.comments = pop_comments(EEG.comments,'','Marked artifact IC components',1);

    end

    EEG = eeg_checkset(EEG);
    %% Clear evenlist and remove events with codelabel '5sec' (created by adjust plugin)

    pop_squeezevents(EEG); % summarize events (ERPLAB)
    % eeg_eventtypes(EEG) % summarize events (EEGLAB)
    EEG.EVENTLIST = [];
    try 
        toIncludeEventsIdx = find(cellfun(@isempty,strfind({EEG.event.type},'5sec')));
        EEG.event = EEG.event(toIncludeEventsIdx);
        pop_squeezevents(EEG); % summarize events (ERPLAB)
        disp('Removed events labeled ''5sec''');
        EEG = eeg_checkset(EEG);
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
    end

    % eeglab redraw
    
    %% Save dataset

    if useICA
        EEG.setname = [EEG.subject '_continuous_ICAweights'];
    else
        EEG.setname = [EEG.subject '_continuous'];
    end
    EEG = pop_saveset(EEG,'filename',EEG.setname,'filepath',fullfile(currentSubjectDirectory,'continuous')); % save EEG set
    if ~isempty(toRemove)
        EEG_emg = pop_saveset(EEG_emg,'filename',[EEG.subject '_emg'],'filepath',fullfile(currentSubjectDirectory,'continuous')); % save EMG set
    end
    disp(['Finished subject ' subject, '!']);
    dlmwrite(['Finished_Subject_' subject '.txt'],'Finished!','delimiter','');

catch ME % if errors when running script, catch them

    disp(['Error with subject ' subject, '!']);
    save(fullfile(pwd, ['Error_MException_' subject '.mat']), 'ME'); % save MException object
    
end % end try catch

clear
close all

end % end  preprocessEEG function
