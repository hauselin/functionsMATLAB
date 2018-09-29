%% Check ICA components and remove blink components one subject at a time.
% Last modified by Hause Lin 15-08-18 13:46 hauselin@gmail.com

clear, clc, close all

%% Set up parameters and read data

subject = '010'; % subject id
cd '/Users/Hause/Dropbox/Working Projects/160201 Framing/Analysis/EEG/Scripts/'
appendEMGdata = true; % add EMG data to EEG data (true or false)
plotICLabels = false;
dataDirectory = 'EEGData'; % directory containing subject/participant data
currentDir = pwd;
parentDir = currentDir(1:end-8); % remove 'Scripts' from path
dataDirectoryAbsPath = fullfile(parentDir, dataDirectory); 
currentSubjectDirectory = fullfile(dataDirectoryAbsPath, subject);
directoryToReadFrom = fullfile(currentSubjectDirectory,'continuous'); % directory with raw data
filesInDirectory = dir(directoryToReadFrom); % files in data directory
fileIdx = find(~cellfun(@isempty, strfind({filesInDirectory.name}, 'ICAweights.set'))); % index of file to read
if length(fileIdx) ~= 1 % if more than one raw data found, error
    error('Check subject id or number of files in data directory!');
end
addpath('/users/hause/dropbox/Apps/MATLAB Toolboxes and Packages/eeglab14_1_2b') % add path on local machine
% Read data with ICA weights and plot ICA components
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % run eeglab
clc
EEG = pop_loadset('filename',filesInDirectory(fileIdx).name,'filepath',filesInDirectory(fileIdx).folder);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG,EEG,0);
disp(['Subject ' EEG.subject])

% Create blink reference channel and save in EEG structure (but don't add it to EEG.data)
% chanSO2Idx = find(strcmpi({EEG.chanlocs.labels}, 'SO1'));
% chanIO2Idx = find(strcmpi({EEG.chanlocs.labels}, 'IO1'));
% EEG.VEOGREFChannel = EEG.data(chanSO2Idx, :) - EEG.data(chanIO2Idx, :); % this channel will be added to data AFTER ICA pruning.

eeglab redraw
close all

% Clear evenlist and remove events with codelabel '5sec' (created by adjust plugin)
EEG.EVENTLIST = [];
try 
    toIncludeEventsIdx = find(cellfun(@isempty,strfind({EEG.event.type},'5sec')));
    EEG.event = EEG.event(toIncludeEventsIdx);
    pop_squeezevents(EEG); % summarize events (ERPLAB)
    disp('Removed events labeled ''5sec''');
    EEG = eeg_checkset(EEG);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
end
% pop_squeezevents(EEG); % summarize events (ERPLAB)
% eeg_eventtypes(EEG) % summarize events (EEGLAB)

% convert triggers from string to numeric
EEG = pop_creabasiceventlist(EEG,'AlphanumericCleaning','on','BoundaryNumeric',{-99},'BoundaryString',{'boundary'});
pop_squeezevents(EEG); % summarize events (ERPLAB)

if isempty(EEG.icaact) % recompute IC time series if it's empy
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    disp('Recomputing IC time series');
end

% plot components
if plotICLabels
    if ~isfield(EEG.etc, 'ic_classification')
        disp('Running iclabel...');
        EEG = iclabel(EEG);
    end
    disp('Done running iclabel...');

    if EEG.nbchan > 40 
        pop_viewprops(EEG,0,1:28,{'freqrange',[2 80]},{},1,'ICLabel'); % plot only 28 components
    else
        pop_viewprops(EEG,0,1:size(EEG.icaact,1),{'freqrange',[2 80]},{},1,'ICLabel'); % plot all components
    end
    print(gcf,'-djpeg','-r100', fullfile(currentSubjectDirectory,'parameters',[EEG.subject '_ICAcomponents.jpg']));
end
clc

% read components removed text file
txtFile = fullfile(currentSubjectDirectory,'parameters',[EEG.subject 'ComponentsRemoved.txt']);
if exist(txtFile)
    fileId = fopen(txtFile,'r');
    componentsToRemove = fscanf(fileId,'%f')';
    fclose(fileId);
    disp('Components have been identified previously!');
    componentsToRemove
end

% read manual bad channel removed mat
matFile = fullfile(currentSubjectDirectory,'parameters',[EEG.subject 'badChannelManualRemove.mat']);
if exist(matFile)
    load(matFile)
    disp('Bad channels were interpolated manually previously!');
    badchannelmanual
end

%% Check IC components metrics

find(EEG.etc.ic_classification.ICLabel.classifications(:, 3) > 0.3 & EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.01) % eye but not brain 
find(EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.0001) % not brain

EEG.etc.ic_classification.ICLabel.classes % iclabel
find(EEG.etc.ic_classification.ICLabel.classifications(:, 2) > 0.95 & EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.015) % muscle
find(EEG.etc.ic_classification.ICLabel.classifications(:, 4) > 0.95 & EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.015) % heart
find(EEG.etc.ic_classification.ICLabel.classifications(:, 6) > 0.95 & EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.015) % channel noise
find(EEG.etc.ic_classification.ICLabel.classifications(:, 5) > 0.95 & EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.015) % line noise
% EEG.adjustresults % adjust results
% EEG.icaquant % icablinkmetrics results
% EEG.etc.ic_classification.ICLabel.classifications

%% Temporarily remove ICA components and compare results

% EEG.componentsRemoved = [1 2 3]; % components to remove
EEG.componentsRemoved = componentsToRemove; % components to remove (from previously saved txt file)
EEG2 = pop_subcomp(EEG,EEG.componentsRemoved,0); % remove components

% plot 
try % close figures if already opened
    close 'EEG raw'
    close 'EEG cleaned (ICA components removed)'
end

% channels to plot to compare pre/post ICA component rejection
toPlot = find(ismember({EEG.chanlocs.labels}, {'Fp1','Fp2','Fpz','FCz','Cz','Fz','Oz','Pz','F7','F8','O1','O2','Oz','IO1','SO1','IO2','IO2','T3','T4','F3','F4','T7','T8','FT7','FT8','SO2','IO2','M1','M2','CPz'}));
scaleRange = 25; % y axis (voltage) range

% plot raw data (all channels)
eegplot(EEG.data(1:EEG.nbchan, :, :),'srate',EEG.srate,'title','EEG raw','eloc_file',EEG.chanlocs,'spacing',scaleRange,'events',EEG.event,'tag','childEEG','position',[10 400 800 700])
% plot cleaned data (all channels)
eegplot(EEG2.data(1:EEG.nbchan, :, :),'srate',EEG2.srate,'title','EEG cleaned (ICA components removed)','eloc_file',EEG2.chanlocs,'spacing',scaleRange,'events',EEG.event,'children',findobj('tag','childEEG'),'position',[800 400 800 700])

% plot raw data (subset of channels)
% eegplot(EEG.data(toPlot, :, :),'srate',EEG.srate,'title','EEG raw','eloc_file',EEG.chanlocs(toPlot),'spacing',scaleRange,'events',EEG.event,'tag','childEEG','position',[10 400 800 700])
% plot cleaned data (subset of channels)
% eegplot(EEG2.data(toPlot, :, :),'srate',EEG2.srate,'title','EEG cleaned (ICA components removed)','eloc_file',EEG2.chanlocs(toPlot),'spacing',scaleRange,'events',EEG.event,'children',findobj('tag','childEEG'),'position',[800 400 800 700])

% plot ICA component time series
% eegplot(EEG.icaact(1:length(EEG.chanlocs), :, :),'srate',EEG.srate,'title','ICA component time series','spacing',50,'events',EEG.event,'tag','childEEG','position',[900 50 800 500])
% close all

%% Remove ICA components permanently

EEG = pop_subcomp(EEG,EEG.componentsRemoved,0); % remove components
disp(['Removed components ' num2str(EEG.componentsRemoved)]);
dlmwrite(fullfile(currentSubjectDirectory,'parameters',[EEG.subject 'ComponentsRemoved.txt']), EEG.componentsRemoved,'delimiter',' ');
EEG.comments = pop_comments(EEG.comments,'','Removed ICA artifact components.',1);
clear EEG2 
close all
% check data
scaleRange = 25;
eegplot(EEG.data(1:EEG.nbchan,:,:),'srate',EEG.srate,'title','EEG cleaned and rereferenced','eloc_file',EEG.chanlocs,'spacing',scaleRange,'events',EEG.event,'position',[10 400 1000 800])

%% Filter and check data quality

% EEG2 = pop_basicfilter(EEG,1:EEG.nbchan,'Boundary','boundary','Cutoff',20,'Design','butter','Filter','lowpass','Order',2,'RemoveDC','off');
% pop_eegplot(EEG2, 1, 1, 1); % plot EEG
% clear EEG2

%% Interpolate bad channels if necessary

badchannelmanual = {}; % manually interpolate bad channel
if ~isempty(badchannelmanual)
   EEG.badchannelmanual = badchannelmanual;
   for chanI = 1:length(badchannelmanual)
       disp(['Interpolating ' badchannelmanual{chanI} '...']);
       EEG = eeg_interp(EEG,find(ismember({EEG.chanlocs.labels},badchannelmanual{chanI}))); % find channel number to interpolate
   end
   EEG.comments = pop_comments(EEG.comments,'','Interpolated bad channels',1);
   disp('Finished interpolating');
   save(fullfile(currentSubjectDirectory,'parameters',[EEG.subject 'badChannelManualRemove.mat']),'badchannelmanual');
else
    disp('No need to interpolate');
end

% add and interpolate channels that have been removed prior to ICA
if isfield(EEG, 'badchannel')
    for chanI=1:length(EEG.badchannel) % add bad channel (as flat channel) back to data
        disp(['Adding ' EEG.badchannel{chanI} '...']);
        EEG.data(end+1,:) = 0; % add empty channel
        if ~isempty(EEG.chanlocs)
            EEG.chanlocs(end+1).labels = EEG.badchannel{chanI}; % add label for channel
        end
    end
    EEG.nbchan = size(EEG.data,1);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
    EEG = pop_chanedit(EEG,'lookup','standard-10-5-cap385.elp'); % add channel locations
    
    % interpolate channels
    chansToInterpolate = (EEG.nbchan-length(EEG.badchannel)+1):EEG.nbchan;
    disp('Interpolating bad channels...');
    % {EEG.chanlocs(chansToInterpolate).labels}
    EEG = eeg_interp(EEG,chansToInterpolate);
    EEG.comments = pop_comments(EEG.comments,'','Added and interpolated bad channels',1);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
    disp('Finished adding and interpolating');
else
    disp('No need to add or interpolate');
end

%% Rereference

% find index for reference channels
ref1Indx = find(strcmpi({EEG.chanlocs.labels}, 'M1'));
ref2Indx = find(strcmpi({EEG.chanlocs.labels}, 'M2'));

% compute and add new channel (average of M1 and M2)
refchan = (EEG.data(ref1Indx,:) + EEG.data(ref2Indx,:)) / 2;

% rereference to refchan
EEG.data = bsxfun(@minus,EEG.data,refchan);
EEG.ref = 'Mastoids average';
EEG.comments = pop_comments(EEG.comments,'','Rereferenced to mastoids averaged',1);
disp('Rereferenced to mastoids averaged');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
EEG = eeg_checkset(EEG);

% final check
scaleRange = 25;
eegplot(EEG.data(1:EEG.nbchan,:,:),'srate',EEG.srate,'title','EEG cleaned and rereferenced','eloc_file',EEG.chanlocs,'spacing',scaleRange,'events',EEG.event,'position',[10 400 1000 800])

%% Read in EMG data, reference them, and add to EEG data

if appendEMGdata
    fileIdx = find(~cellfun(@isempty, strfind({filesInDirectory.name}, '_emg.set'))); % index of file to read
    if length(fileIdx) > 1 % if more than one raw data found, error
        error('Too many matching file names in data directory!');
    elseif length(fileIdx) == 0
        error('No emg .set files found!');
    elseif length(fileIdx) == 1
        disp('Reading in and processing EMG data...');
        % Read data with ICA weights and plot ICA components
        EEG = pop_loadset('filename',filesInDirectory(fileIdx).name,'filepath',filesInDirectory(fileIdx).folder);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG);

        if (ALLEEG(1).pnts ~= ALLEEG(2).pnts)
            error('Number of data points in EEG and EMG data do not match!');
        end
        
        {EEG.chanlocs.labels} % all emg channel labels

        % 28 to 500 Hz bandpass filter
        EEG = pop_basicfilter(EEG,1:EEG.nbchan,'Boundary','boundary','Cutoff',[28 500],'Design','butter','Filter','bandpass','Order',4,'RemoveDC','on'); 
        % notch filter to remove line noise
        EEG = pop_basicfilter(EEG,1:EEG.nbchan,'Boundary','boundary','Cutoff',60,'Design','notch','Filter','PMnotch','Order',180); 
        disp('Filtered EMG data.');

        % compute rereferenced corrugator channel
        EMG.data = EEG.data(strcmpi('CORRins',{EEG.chanlocs.labels}),:) - EEG.data(strcmpi('CORRout',{EEG.chanlocs.labels}),:);
        EMG.chan = {'Corr'};

        % compute rereferenced zygomat channel
        EMG.data(end+1,:) = EEG.data(strcmpi('ZYGlow',{EEG.chanlocs.labels}),:) - EEG.data(strcmpi('ZYGup',{EEG.chanlocs.labels}),:);
        EMG.chan(end+1) = {'Zygo'};
        disp('Rereferenced EMG data.');

        % calculate weights for convolving (moving average) (filter data)
        windowMs = 20; % moving average window in ms
        windowPoints = EEG.srate / 1000 * windowMs; % moving average in data points
        windowPoints = round(windowPoints);
        if mod(windowPoints, 2) == 0 % if window size is even, make it odd by adding 1 point
            windowPoints = windowPoints + 1;
        end
        weights = ones(windowPoints, 1) / windowPoints; % weights for convolution
        
        % rectify (take absolute value of reref channel) and convolve
        for i=1:length(EMG.chan)
            EMG.data(i,:) = abs(EEG.data(i,:));
            EMG.data(i,:) = conv(EMG.data(i,:), weights,'same'); % convolve/filter data using moving average
        end
        disp('Convolved (moving average) EMG data.');
        disp('Rectified (take absolute value) EMG data.');

        EEG = eeg_retrieve(ALLEEG,1); CURRENTSET = 1; % retrieve EEG dataset
        EEG = eeg_checkset(EEG);

        % add rereferenced EMG channels back to EEG data
        for i=1:length(EMG.chan) 
            EEG.data(end+1,:) = EMG.data(i,:);
            EEG.nbchan = size(EEG.data,1);
            EEG.chanlocs(end+1).labels = EMG.chan{i};
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
        end
        EEG.comments = pop_comments(EEG.comments,'','Added filtered (28 to 500 Hz) EMG channels.',1);
        disp('Added filtered (28 to 500 Hz) EMG channels to EEG data.');
    end
    EEG = eeg_checkset(EEG);
    EEG = pop_chanedit(EEG,'lookup','standard-10-5-cap385.elp');

end
close all

%% Reorder channels

% eegplot(EEG.data(1:EEG.nbchan,:,:),'srate',EEG.srate,'title','EEG cleaned and rereferenced','eloc_file',EEG.chanlocs,'spacing',scaleRange,'events',EEG.event,'position',[10 400 1000 800])
newChanOrder = {'Fpz' 'Fz' 'FCz' 'Cz' 'CPz' 'Pz' 'Oz' 'M1' 'M2' 'SO2' 'IO2' 'Corr' 'Zygo'};
EEG = reorderchans(EEG,newChanOrder,true);
print(gcf,'-djpeg','-r100', fullfile(currentSubjectDirectory,'parameters',[EEG.subject 'chanLocs.jpg']));
EEG = eeg_checkset(EEG);
disp('Reordered channels.');

% final check
scaleRange = 30;
eegplot(EEG.data(1:EEG.nbchan,:,:),'srate',EEG.srate,'title','EEG cleaned, rereferenced, sorted','eloc_file',EEG.chanlocs,'spacing',scaleRange,'events',EEG.event,'position',[10 400 1000 800])

%% save data

EEG.setname = [EEG.subject '_ICAPruned_Reref'];
EEG = pop_saveset(EEG,'filename',EEG.setname,'filepath',EEG.filepath);
disp(['Saved ' fullfile(EEG.filepath,EEG.setname)]);

clear, clc
close all