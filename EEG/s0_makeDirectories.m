% Copies raw data to a new directory with the following structure (example):
% 
% |-- EEGAnalysis
%   |-- Binlister
%   |-- EEGData
%       |-- 001
%       |-- 002
%           |-- continuous
%           |-- parameters
%           |-- raw
%               |-- Stroop002.cnt
%               |-- Stroop002.trg
%       |-- 003
%       |-- 004
%   |-- Figures
%   |-- Messages
%   |-- README.md
%   |-- Scripts
%
% USAGE 
% s0_makeDirectories({'0001','0002'}) % no directories provided (pop-up GUI to first select input then output directories)
% s0_makeDirectories('001','/Users/Dropbox/Data/MyRawData/','/Users/Experiments/') % one subject
% s0_makeDirectories({'0001','0002','0003','0013'},'/Users/DataFromEEG/','/Users/Experiments/') % multiple subjects at once
% s0_makeDirectories({'0001','0002','0003','0013'},'Raw','/Users/Experiments/') % raw data in subdirectory 'Raw'
%
% INPUTS
% subjectIDs: subject IDs STRING to match in your raw data filenames (can also be a CELL array containing multiple subject IDs)
%             partial matching works too: 'EEGStudy001.edf' and 'EEGStudy020.edf' can be matched using just {'001','020'}
% rawDataDirectory: raw data directory STRING (path to raw data)
% outputDirectory: output directory STRING

% Written in MATLAB R2016b (also tested in 2017b)
% Last modified by Hause Lin 04-08-18 08:37 hauselin@gmail.com

function s0_makeDirectories(subjectIDs,rawDataDirectory,outputDirectory)

%% set up new data directory structure

mainDirectory = 'EEGAnalysis'; % new directories to make
directoriesToMake = {'Scripts' 'Messages' 'EEGData' 'Binlister' 'Results' 'Figures'}; % directories in main directory
newSubjectDirectoriesToMake = {'continuous','raw','parameters'}; % subdirectories to be created within each participant's folder

%% manually generate abd overwrite subjectIDs to match if necessary... (overwrites input subjectIDs)

% subjectIDs = {'1.edf' '12.edf'};
% create cleaned output subject numbers
% subjectIDsOutput = cellfun(@str2num,(strrep(subjectIDs,'.edf','')));
% subjectIDsOutput = strsplit(num2str(subjectIDsOutput, '%03d '), ' '); % convert numeric to string in cell to match file names and move files accordingly

% subjectIDs = 1:100;
% subjectIDs = strsplit(num2str(subjectIDs, '%03d '), ' '); % convert numeric to string in cell to match file names and move files accordingly

%% Verifying input/output before running

if isa(subjectIDs,'char') % if input is character (e.g., '001'), convert to cell array (e.g., {'001'})
   subjectIDs = {subjectIDs};
end
subjectIDsOutput = subjectIDs; 

clc;
if nargin == 1
    disp('Where to read files from?');
    rawDataDirectory = uigetdir(pwd,'Raw file location?');
    disp(['Reading raw data from ' rawDataDirectory]);
    disp('Where to copy files to?');
    outputDirectory = uigetdir(pwd,'Output location?');
    disp(['Copying to this directory: ' fullfile(outputDirectory,mainDirectory)]); 
else
    disp(['Reading raw data from ' rawDataDirectory]);
    disp(['Copying to this directory: ' fullfile(outputDirectory,mainDirectory)]); 
end
input(['If all correct, press enter. If wrong, press Ctrl-C to stop...']);

disp(['Trying to match and move ' num2str(length(subjectIDs)) ' subjects']);
disp(subjectIDs);
disp(['New subject output directories']);
disp(fullfile('EEGData',subjectIDsOutput));
disp('Creating these directories within each subject''s directory:');
disp(newSubjectDirectoriesToMake)
input(['If all correct, press enter. If wrong, press Ctrl-C to stop...']);

for dirI = 1:length(directoriesToMake)
    if ~exist(fullfile(outputDirectory,mainDirectory,directoriesToMake{dirI}))
        mkdir(fullfile(outputDirectory,mainDirectory,directoriesToMake{dirI}));
    end
end

if ~exist(fullfile(outputDirectory,mainDirectory,'README.md')) % create README.md file if doesn't exist
    dlmwrite(fullfile(outputDirectory,mainDirectory,'README.md'),'# Information about analysis','delimiter','');
end

%% loop through each subject, create directories, subdirectories, and move files into subdirectories

clc
for subjI = 1:length(subjectIDs)
    
    disp('---------------------------------------------------------');
    currentSubject = subjectIDs{subjI}; % get current subject id
    disp(['Looking for ' currentSubject]);
    
    % find matching files to move in top directory (rawDataDirectory)
    filesInRawDataDirectory = dir(rawDataDirectory);
    filesToMove = find(~cellfun(@isempty, strfind({filesInRawDataDirectory.name},currentSubject))); % indices of files to read (partial match)
    % filesToMove = find(ismember({filesInRawDataDirectory.name},currentSubject)); % indices of files to read (exact match)
    
    if ~isempty(filesToMove) % if matching files found
        
        % create necessary directories for subject inside EEGData/subjectID directory
        subjDirectory = fullfile(outputDirectory,mainDirectory,'EEGData',subjectIDsOutput{subjI});
        for dirI = 1:length(newSubjectDirectoriesToMake)
            if ~exist(fullfile(subjDirectory,newSubjectDirectoriesToMake{dirI})) % if don't exist, create
                mkdir(fullfile(subjDirectory,newSubjectDirectoriesToMake{dirI}));
            end
        end

        moveFiles = {filesInRawDataDirectory(filesToMove).name}; % cell array of files to move
        for fileI = 1:length(moveFiles) % move files into 'EEGData/subjectID/raw' folder
            copyfile(fullfile(rawDataDirectory,moveFiles{fileI}),fullfile(subjDirectory,'raw')); % copy files
            disp(['Copied ' moveFiles{fileI}]);
        end
        disp(['Finished subject ' num2str(subjectIDs{subjI})]);
        disp('---------------------------------------------------------');
        
    else
        
        disp(['Failed to find any files for subject ' num2str(subjectIDs{subjI}) '!']);
        disp('---------------------------------------------------------');
        
    end

end

end