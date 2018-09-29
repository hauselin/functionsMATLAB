function preprocessEMG(subject)
% Preprocesses EMG data.
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

eegdataformat = '.cnt'; % raw data file format
dataDirectory = 'EEGData'; % directory containing subject/participant data

% trim/remove long sections of data without events (e.g., between-block breaks)
minTimeThreshold = 5000; % minimum time (s) between two events
bufferTime = 4500; % add buffer time (s) (pre/post event)

% emg channel (to be removed from this preprocessing)
emgChan = {'CorsOut','CorsIns','CORRins','CORRout','ZYGup','ZYGlow','COORins','COORout'}; % name of emg channels

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

    % remove emg channels
    if ~isempty(emgChan)
        toRemove = [];
        for emgChanI = emgChan
            toRemove = [toRemove find(strcmpi(emgChanI{1}, allChannels))]; % find index of emg channels
        end
        toRemove = allChannels(toRemove); % channel labels to remove
        EEG = pop_select(EEG,'channel',toRemove); % save EMG channels
    end

    %% Trim data (ERPLAB function)

    % trim EEG data
    EEG = pop_erplabDeleteTimeSegments(EEG,'displayEEG',0,'endEventcodeBufferMS',bufferTime,'ignoreUseType','ignore','startEventcodeBufferMS',bufferTime,'timeThresholdMS',minTimeThreshold);
    EEG.comments = pop_comments(EEG.comments,'','Trimmed data by deleting eventless time segments.',1);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,1); % modify ALLEEG structure; check consistency; store in dataset1
    
    %% Checking events

    % pop_squeezevents(EEG); % summarize events (ERPLAB)
    % eeg_eventtypes(EEG) % summarize events (EEGLAB)
    % convert triggers from string to numeric
    EEG = pop_creabasiceventlist(EEG,'AlphanumericCleaning','on','BoundaryNumeric',{-99},'BoundaryString',{'boundary'});
    EEG = eeg_checkset(EEG);
    %% Clear evenlist and remove events with codelabel '5sec' (created by adjust plugin)

    pop_squeezevents(EEG); % summarize events (ERPLAB)
    % eeg_eventtypes(EEG) % summarize events (EEGLAB)
    % eeglab redraw
    
    %% Save dataset

    EEG.setname = [EEG.subject '_emg'];
    EEG = pop_saveset(EEG,'filename',EEG.setname,'filepath',fullfile(currentSubjectDirectory,'continuous')); % save EEG set

    disp(['Finished subject ' subject, '!']);
    dlmwrite(['Finished_Subject_' subject '.txt'],'Finished!','delimiter','');

catch ME % if errors when runnig script, catch them

    disp(['Error with subject ' subject, '!']);
    save(fullfile(pwd, ['Error_MException_' subject '.mat']), 'ME'); % save MException object

end % end try catch

clear
close all

end % end  preprocessEEG function
