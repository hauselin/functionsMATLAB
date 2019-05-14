%% Check ICA components and remove blink components one subject at a time.
% Last modified by Hause Lin 15-08-18 13:46 hauselin@gmail.com

eeglab
clear, clc, close all
cd '/Users/hause/Dropbox/Working Projects/Akina Physical Effort/Zaibas RewP study/Scripts/'

%% Set up paths and read data

subject = '3'; % subject id

PATHS = struct();
PATHS.datadir = '/Users/hause/Dropbox/Working Projects/Akina Physical Effort/Zaibas RewP study/DataPreprocessed/';
PATHS.outputdir = '/Users/hause/Dropbox/Working Projects/Akina Physical Effort/Zaibas RewP study/DataICA_cleaned/';

cfg = struct();
cfg.plotICLabels = true;

load(fullfile(PATHS.datadir,[subject '_continuous_icaed.mat']))
clc; disp(['Subject ' EEG.subject]);

% ensure triggers are numeric
EEG = pop_creabasiceventlist(EEG,'AlphanumericCleaning','on','BoundaryNumeric',{-99},'BoundaryString',{'boundary'});
% pop_squeezevents(EEG); % summarize events (ERPLAB)
% eeg_eventtypes(EEG) % summarize events (EEGLAB)

if isempty(EEG.icaact) % recompute IC time series if it's empy
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
end

% plot components
if cfg.plotICLabels
    if ~isfield(EEG.etc, 'ic_classification')
        EEG = iclabel(EEG);
    end
    if EEG.nbchan > 35
        pop_viewprops(EEG,0,1:35,{'freqrange',[2 70]},{},1,'ICLabel'); % plot only 28 components
    else
        pop_viewprops(EEG,0,1:size(EEG.icaact,1),{'freqrange',[2 70]},{},1,'ICLabel'); % plot all components
    end
end
clc

%% Check IC components metrics

find(EEG.etc.ic_classification.ICLabel.classifications(:, 3) > 0.3 & EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.01) % eye but not brain 
find(EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.0001) % not brain
EEG.etc.ic_classification.ICLabel.classes % iclabel
find(EEG.etc.ic_classification.ICLabel.classifications(:, 2) > 0.95 & EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.015) % muscle
find(EEG.etc.ic_classification.ICLabel.classifications(:, 4) > 0.95 & EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.015) % heart
find(EEG.etc.ic_classification.ICLabel.classifications(:, 6) > 0.95 & EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.015) % channel noise
find(EEG.etc.ic_classification.ICLabel.classifications(:, 5) > 0.95 & EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.015) % line noise
% EEG.icaquant % icablinkmetrics results

%% Temporarily remove ICA components; plot and compare results

EEG.componentsRemoved = [1 4 6 9]; % components to remove
EEG2 = pop_subcomp(EEG,EEG.componentsRemoved,0); % remove components

% plot 
try % close figures if already opened
    close 'EEG raw'
    close 'EEG cleaned (ICA components removed)'
end

% channels to plot to compare pre/post ICA component rejection
toPlot = find(ismember({EEG.chanlocs.labels}, {'Fp1','Fp2','Fpz','FCz','Cz','Fz','Oz','Pz','F7','F8','O1','O2','Oz','IO1','SO1','IO2','IO2','T3','T4','F3','F4','T7','T8','FT7','FT8','SO2','IO2','M1','M2','CPz'}));
scaleRange = 35; % y axis (voltage) range

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

%% Permanently remove ICA components

EEG = pop_subcomp(EEG,EEG.componentsRemoved,0); % remove components
disp(['Removed components ' num2str(EEG.componentsRemoved)]);
EEG.comments = pop_comments(EEG.comments,'','Removed ICA artifact components.',1);
clear EEG2 
close all
% check data
scaleRange = 25;
eegplot(EEG.data(1:EEG.nbchan,:,:),'srate',EEG.srate,'title','EEG cleaned and rereferenced','eloc_file',EEG.chanlocs,'spacing',scaleRange,'events',EEG.event,'position',[10 400 1000 800])

%% Check data quality of lowpass filtered data

EEG2 = pop_basicfilter(EEG,1:EEG.nbchan,'Boundary','boundary','Cutoff',30,'Design','butter','Filter','lowpass','Order',2,'RemoveDC','off');
pop_eegplot(EEG2, 1, 1, 1); % plot EEG
clear EEG2

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
    EEG = pop_chanedit(EEG,'lookup','standard-10-5-cap385.elp'); % add channel locations
    
    % interpolate channels
    chansToInterpolate = (EEG.nbchan-length(EEG.badchannel)+1):EEG.nbchan;
    disp('Interpolating bad channels...');
    % {EEG.chanlocs(chansToInterpolate).labels}
    EEG = eeg_interp(EEG,chansToInterpolate);
    EEG.comments = pop_comments(EEG.comments,'','Added and interpolated bad channels',1);
    disp('Finished adding and interpolating');
else
    disp('No need to add or interpolate');
end

%% Rereference

% find index for reference channels
ref1Indx = find(strcmpi({EEG.chanlocs.labels}, 'M1'));
ref2Indx = find(strcmpi({EEG.chanlocs.labels}, 'M2'));
refchan = (EEG.data(ref1Indx,:) + EEG.data(ref2Indx,:)) / 2; % compute and add new channel (average of M1 and M2)

% rereference to refchan
EEG.data = bsxfun(@minus,EEG.data,refchan);
EEG.ref = 'Mastoids average';
EEG.comments = pop_comments(EEG.comments,'','Rereferenced to mastoids averaged',1);
disp('Rereferenced to mastoids averaged');
EEG = eeg_checkset(EEG);

% final check
scaleRange = 25;
eegplot(EEG.data(1:EEG.nbchan,:,:),'srate',EEG.srate,'title','EEG cleaned and rereferenced','eloc_file',EEG.chanlocs,'spacing',scaleRange,'events',EEG.event,'position',[10 400 1000 800])

%% Reorder channels

newChanOrder = {'Fp1','Fpz','Fp2','F8','F4','Fz','F3','F7','FC3','FCz','FC4','C4','Cz','C3','FT7','T3','FT8','T4','CP4','CPz','CP3','P3','Pz','P4','TP8','T6','TP7','T5','O1','Oz','O2','M1','M2','SO2','IO2','LO1','LO2'};
EEG = reorderchans(EEG,newChanOrder,true);
EEG = eeg_checkset(EEG);
disp('Reordered channels.');

% final check
scaleRange = 30;
eegplot(EEG.data(1:EEG.nbchan,:,:),'srate',EEG.srate,'title','EEG cleaned, rereferenced, sorted','eloc_file',EEG.chanlocs,'spacing',scaleRange,'events',EEG.event,'position',[10 400 1000 800])

%% save data

if ~exist(PATHS.outputdir)
    mkdir(PATHS.outputdir);
end
save(fullfile(PATHS.outputdir,[subject '_continuous_icacleaned.mat']),'EEG');

clear, clc
close all