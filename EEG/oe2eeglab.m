function EEG = oe2eeglab(pathname,filebasename)
%
% EEG = oe2eeglab(pathname,basefilename)
%
%   Loads open-ephys-format data into EEGLAB structure.
%
%   Inputs:
%     pathname: the path to the folder where the files are.
%     filename: any filename up to '.continuous'
% 
%   (if no inputs are given, a dialog box will open to select a file)
%
%   Outputs:
%
%     EEG: EEG structure (see eeglab)
%          *note* data are downsampled to 1 kHz
% 
% 
% hacked from open-ephys original code by mikexcohen@gmail.com

%
%   DISCLAIMER:
%
%   Both the Open Ephys data format and this m-file are works in progress.
%   There's no guarantee that they will preserve the integrity of your
%   data. They will both be updated rather frequently, so try to use the
%   most recent version of this file, if possible.
%
%

%
%     ------------------------------------------------------------------
%
%     Copyright (C) 2014 Open Ephys
%
%     ------------------------------------------------------------------
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     <http://www.gnu.org/licenses/>.
%


%% initial stuff

% no inputs, select a file
if nargin==0
    [filebasename,pathname] = uigetfile('*.continuous','Select any channel file');
end


% housekeeping /\
if ~isequal(pathname(end),filesep), pathname = [pathname filesep]; end

disp(' ')
disp([ 'Importing data from ' pathname filebasename '...' ])
disp(' ')

% get list of data channels to import
underscore = strfind(filebasename,'CH');
allfiles = dir([ pathname filebasename(1:underscore-1) 'CH*.continuous' ]);
ADCfiles = dir([ pathname filebasename(1:underscore(end)) '*ADC*.continuous' ]);


% remove unwanted children
tokeep = true(size(allfiles));
% if ~isequal(filebasename(end-1),'_')
%     for i=1:length(allfiles)
%         if isequal(allfiles(i).name(end-12),'_')
%             tokeep(i) = 0;
%         end
%     end
% end

% change format
allfiles = {allfiles(tokeep).name};

% resort according to number instead of string
fnum = zeros(size(allfiles));
for i=1:length(allfiles)
    u = strfind(allfiles{i},'CH');
    u = u(1)+2;
    d = [strfind(allfiles{i}(u:end),'.') strfind(allfiles{i}(u:end),'_') ];
    fnum(i) = sscanf(allfiles{i}(u:min(d)+u-2),'%g');
end
[~,ford] = sort(fnum);
allfiles = allfiles(ford);

% prepare empty EEG structure
EEG = struct('setname',[],'filename',[],'filepath',[],'pnts',[],'nbchan',[],'trials',[],'srate',[],'xmin',[],'xmax',[],'data',[],'icawinv',[],'icasphere',[],'icaweights',[],'icaact',[],'event',[],'epoch',[],'chanlocs',[],'chaninfo',[],'comments',[],'ref',[],'saved',[]);

%% loop over channels and import data, downsampling to 1 kHz

for filei = 1:length(allfiles)
    
    %% import from raw
    
    % open file
    fid = fopen([pathname allfiles{filei}]);
    fseek(fid,0,'eof');
    filesize = ftell(fid);
    
    % read header
    NUM_HEADER_BYTES = 1024;
    fseek(fid,0,'bof');
    info = getHeader( fread(fid, NUM_HEADER_BYTES, 'char*1') );
    dblock = struct('Repeat',{1 1 1 1024 10},'Types', {'int64' 'uint16' 'uint16' 'int16' 'uint8'},'Str',{'ts' 'nsamples' 'recNum' 'data' 'recordMarker'});
    
    blockBytes = str2double(regexp({dblock.Types},'\d{1,2}$','match', 'once')) ./8 .* cell2mat({dblock.Repeat});
    numIdx = floor((filesize - NUM_HEADER_BYTES)/sum(blockBytes));
    
    % now read data
    info.ts = segRead('ts');
    info.nsamples = segRead('nsamples');
    data = segRead('data','b').*info.header.bitVolts; % read in data
    timestamps = nan(size(data));
    current_sample = 0;
    for record = 1:length(info.ts)
        timestamps(current_sample+1:current_sample+info.nsamples(record)) = info.ts(record):info.ts(record)+info.nsamples(record)-1;
        current_sample = current_sample + info.nsamples(record);
    end
    fclose(fid);
    
    %% downsample and put in EEG structure
    
    if filei==1
        EEG.data = zeros(length(allfiles),ceil(current_sample/30));
    end
    fprintf('.');
    
    % anti-aliasing filter
    ford = 30*(info.header.sampleRate/500);
    fkern = fir1(ford,500/(info.header.sampleRate/2));
    
    % zero-phase-shift filter with reflection
    dataR = [data(ford:-1:1); data; data(end:-1:end-ford+1)]; % reflect
    data = filter(fkern,1,dataR);           % forward filter
    data = filter(fkern,1,data(end:-1:1));  % reverse filter
    data = data(end:-1:1);                  % reverse again for 0phase
    data = data(ford+1:end-ford);           % chop off reflected parts
    
    % put downsampled data into EEG structure
    EEG.data(filei,:) = data(1:30:end);
    EEG.chanlocs(filei).labels = num2str(filei);
    
end

EEG.srate = info.header.sampleRate/30;
EEG.times = timestamps(1:30:end) / info.header.sampleRate;

%% Import analogue channel if it exists

if ~isempty(ADCfiles)
    ADCfiles = {ADCfiles.name};
    
    for filei = 1 %:length(ADCfiles) Change this if there is more than Analogue channel
        
        %% import from raw
        
        % open file
        fid = fopen([pathname ADCfiles{filei}]);
        fseek(fid,0,'eof');
        filesize = ftell(fid);
        
        % read header
        NUM_HEADER_BYTES = 1024;
        fseek(fid,0,'bof');
        info = getHeader( fread(fid, NUM_HEADER_BYTES, 'char*1') );
        dblock = struct('Repeat',{1 1 1 1024 10},'Types', {'int64' 'uint16' 'uint16' 'int16' 'uint8'},'Str',{'ts' 'nsamples' 'recNum' 'data' 'recordMarker'});
        
        blockBytes = str2double(regexp({dblock.Types},'\d{1,2}$','match', 'once')) ./8 .* cell2mat({dblock.Repeat});
        numIdx = floor((filesize - NUM_HEADER_BYTES)/sum(blockBytes));
        
        % now read data
        info.ts = segRead('ts');
        info.nsamples = segRead('nsamples');
        data = segRead('data','b').*info.header.bitVolts; % read in data
        timestamps = nan(size(data));
        current_sample = 0;
        for record = 1:length(info.ts)
            timestamps(current_sample+1:current_sample+info.nsamples(record)) = info.ts(record):info.ts(record)+info.nsamples(record)-1;
            current_sample = current_sample + info.nsamples(record);
        end
        fclose(fid);
        
        %% downsample and put in EEG structure
        
        EEG.data(end+filei,:) = data(1:30:end);
        EEG.chanlocs(end+filei).labels = ['ADC_' num2str(filei)];
        
    end
end
%% import markers

% find out if there's a _* in the filename
% und = strfind(filebasename,'_');
% if numel(und)==1
    efile = dir([ pathname 'all_channels.events' ]);
% else
%     efile = dir([ pathname 'all_channels_' filebasename(end) '.events' ]);
% end


fid = fopen([ pathname efile.name ]);
fseek(fid,0,'eof');
filesize = ftell(fid);

NUM_HEADER_BYTES = 1024;
fseek(fid,0,'bof');
hdr = fread(fid, NUM_HEADER_BYTES, 'char*1');
info = getHeader(hdr);

dblock = struct('Repeat',{1 1 1 1 1 1 1},'Types', {'int64' 'uint16' 'uint8' 'uint8' 'uint8' 'uint8' 'uint16'},'Str',{'timestamps' 'sampleNum' 'eventType' 'nodeId' 'eventId' 'data' 'recNum'});
blockBytes = str2double(regexp({dblock.Types},'\d{1,2}$','match', 'once')) ./8 .* cell2mat({dblock.Repeat});
numIdx = floor((filesize - NUM_HEADER_BYTES)/sum(blockBytes));
timestamps = segRead('timestamps')./info.header.sampleRate;
info.sampleNum = segRead('sampleNum');
info.eventType = segRead('eventType');
info.nodeId = segRead('nodeId');
info.eventId = segRead('eventId');
data = segRead('data');

%%

data(info.eventType==5) = [];
timestamps(info.eventType==5) = [];
info.eventId(info.eventType==5) = [];
info.eventType(info.eventType==5) = [];

if ~isempty(info.eventId)
    
    event_num = 1;
    event_Id = info.eventId(1);
    event_count = zeros(1,length(data));
    
    for j = 1:length(data)
        if info.eventId(j) == event_Id
            event_count(j) = event_num;
        else
            event_num = event_num + 1;
            event_Id = info.eventId(j);
            event_count(j) = event_num;
        end
    end
    
    for j = 1:event_num
        codes(j) = sum(2.^(data(event_count==j)));
    end
    
    for j = 2:2:event_num
        if codes(j) == codes(j-1)
            events.marker(j./2) = codes(j);
            events.time_on(j./2) = mean(timestamps(event_count==j-1));
            events.time_off(j./2) = mean(timestamps(event_count==j));
        else
            events.time_on(j./2) = NaN;
            events.time_off(j./2) = NaN;
        end
    end
    
    % Check for some weird bug - I think due to dropped samples missing
    % some event data
   
    oddEvents = find(events.marker == 0);
    events.marker(oddEvents) = [];
    events.time_on(oddEvents) = [];
    events.time_off(oddEvents) = [];
    
%% convert to eeglab format - times are given in points...
    
    time_on  = dsearchn(EEG.times,events.time_on'./info.header.sampleRate)';
    time_off = dsearchn(EEG.times,events.time_off'./info.header.sampleRate)';
    
    % time_on  = (events.time_on ./info.header.sampleRate) - EEG.times(1);
    % time_off = (events.time_off ./info.header.sampleRate) - EEG.times(1);
    EEG.event = struct('type',num2cell(events.marker),'latency', num2cell(time_on),'duration',num2cell(time_off - time_on));
    
end


%% other EEG fields

EEG.times = EEG.times-EEG.times(1);

EEG.setname    = 'Raw data';
EEG.filepath   = pathname;
EEG.pnts       = size(EEG.data,2);
EEG.nbchan     = size(EEG.data,1);
EEG.trials     = 1;
EEG.xmin       = EEG.times(1);
EEG.xmax       = EEG.times(end);

fprintf('\n');
%% extra functions

function seg = segRead(segName, mf)
    if nargin == 1, mf = 'l'; end
    segNum = find(strcmp({dblock.Str},segName));
    fseek(fid, sum(blockBytes(1:segNum-1))+NUM_HEADER_BYTES, 'bof'); 
    seg = fread(fid, numIdx*dblock(segNum).Repeat, sprintf('%d*%s', dblock(segNum).Repeat,dblock(segNum).Types), sum(blockBytes) - blockBytes(segNum), mf);
end
end
function info = getHeader(hdr)
    eval(char(hdr'));
    info.header = header;
end

