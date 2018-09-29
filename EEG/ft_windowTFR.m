function [output] = ft_windowTFR(cfg,data)
% Determine time-frequency window
% 
% cfg.xtoi: time in seconds (not necessary if input is fieldtrip structure)
% cfg.yfoi: frequencies (not necessary if input is fieldtrip structure)
% cfg.plotstat: 'minmax' or 'min' or 'max' or 'yes' or 'no' (default 'no'; else requires ft_contourTFR function)
% cfg.subplot: subplot number (default [1 1 1]);
% cfg.ylog = 'yes' (default) or 'no'
% cfg.chan: which channel to plot (if input is ft structure: default is FCz; if not, just plot channel 1); 
% cfg.trials: which trial of the 'rpt' dimension to plot (if input is ft structure; default trial 1)
% cfg.window1_xlim: time range to find min, max, and mean activity in window1 (big window) (default all time points)
% cfg.window1_ylim: frequency range to find min, max, and mean activity in window1 (big window) (default all frequencies)
% cfg.window2_xrange: width of window around min or max activity (default: entire data)
% cfg.window2_yrange: height of window around min or max activity (default: entire data)
% cfg.win2insidewin1: whether window2 must remain inside window1 (default: yes)
% cfg.saveas: specify where to save and file name
%
% Written in MATLAB R2017b
% Last modified by Hause Lin 26-08-18 15:41 hauselin@gmail.com

%% check input data structure and set default values

if ~isstruct(data) && ndims(data) == 2 && ismatrix(data) % if 2D matrix, then just plot the data (provide time and freq: xtoi and yfoi in cfg)
    if ~isfield(cfg,'chan') cfg.chan = 'Channel'; end
    if ~isfield(cfg,'xtoi') error('Specify cfg.xtoi'); end
    if ~isfield(cfg,'yfoi') error('Specify cfg.yfoi'); end
    dataMat = data; % freq_time matrix
elseif isstruct(data) && ndims(data.powspctrm) == 3 % multi-channel data
    if ~isfield(cfg,'chan') cfg.chan = 'FCz'; end
    cfg.chanIndex = find(strcmpi(cfg.chan,data.label));
    if isempty(cfg.chanIndex) 
        cfg.chanIndex = 1; 
        cfg.chan = data.label{cfg.chanIndex};
    end
    cfg.xtoi = data.time;
    cfg.yfoi = data.freq;
    dataMat = squeeze(data.powspctrm(cfg.chanIndex,:,:)); % freq_time matrix
elseif isstruct(data) && ndims(data.powspctrm) == 4 % multi-trial and channel data
    if ~isfield(cfg,'chan') cfg.chan = 'FCz'; end
    cfg.chanIndex = find(strcmpi(cfg.chan,data.label));
    if isempty(cfg.chanIndex)
        cfg.chanIndex = 1; 
        cfg.chan = data.label{cfg.chanIndex};
    end
    cfg.xtoi = data.time;
    cfg.yfoi = data.freq;
    if ~isfield(cfg,'trials') cfg.trials = 1; end % default select only first trial
    dataMat = squeeze(data.powspctrm(cfg.trials,cfg.chanIndex,:,:)); % freq_time matrix
else
    error('Check input data!')
end

% whether to window2 should always be inside window1 (default yes)
if ~isfield(cfg,'win2insidewin1') cfg.win2insidewin1 = 'yes'; end


%% determine window to find min/max value is all data points in 2D matrix

% default is to search within the entire 2D matrix (no windowing)
if ~isfield(cfg,'window1_xlim') 
    cfg.window1_xlim = [cfg.xtoi(1) cfg.xtoi(end)];
    cfg.windowdata = 'no';
else
    cfg.windowdata = 'yes'; % whether to search for min within a window of current 2D mat
end 
if ~isfield(cfg,'window1_ylim') 
    cfg.window1_ylim = [cfg.yfoi(1) cfg.yfoi(end)]; 
    cfg.windowdata = 'no';
else
    cfg.windowdata = 'yes'; % whether to search for min within a window of current 2D mat
end 

% select data window
if strcmpi(cfg.windowdata,'yes')
    dataMatWindow = dataMat;
    % x: convert values outside window to NaN
    cfg.window1_xlimidx = cfg.xtoi >= cfg.window1_xlim(1) & cfg.xtoi <= cfg.window1_xlim(2); % indices inside window
    dataMatWindow(:,~cfg.window1_xlimidx) = NaN;
    % y: convert values outside window to NaN
    cfg.window1_ylimidx = cfg.yfoi >= cfg.window1_ylim(1) & cfg.yfoi <= cfg.window1_ylim(2); % indices inside window
    dataMatWindow(~cfg.window1_ylimidx,:) = NaN;
else
    dataMatWindow = dataMat; % entire dataset
end
% contourf(dataMatWindow,'edgecolor','none')

%% determine values in window 1 (initial window)

output = [];

output.window1.desc = 'x: time, y = freq';
output.window1.xlim = cfg.window1_xlim;
output.window1.ylim = cfg.window1_ylim;
output.window1.rectangle = 'xleft_ybottom_width_yheight';
output.window1.rectangleposition = [cfg.window1_xlim(1) cfg.window1_ylim(1) cfg.window1_xlim(2)-cfg.window1_xlim(1) cfg.window1_ylim(2)-cfg.window1_ylim(1)];

% compute mean value in inital search window (window1)
output.window1.mean = nanmean(dataMatWindow(:));

% find max value in window1
[maxval maxvalidx] = nanmax(dataMatWindow(:)); % find max value and index
[idx1 idx2] = ind2sub(size(dataMatWindow),maxvalidx); % convert vector index to matrix index
output.window1.max = maxval; % max value
output.window1.maxidx = [idx1 idx2];
output.window1.maxidx = [idx1 idx2];
output.window1.maxtoi = [cfg.xtoi(idx2)];
output.window1.maxfoi = [cfg.yfoi(idx1)];

% find min value in window1
[minval minvalidx] = nanmin(dataMatWindow(:));
[idx1 idx2] = ind2sub(size(dataMatWindow),minvalidx);
output.window1.min = minval;
output.window1.minidx = [idx1 idx2];
output.window1.mintoi = [cfg.xtoi(idx2)];
output.window1.minfoi = [cfg.yfoi(idx1)];

% dataMatWindow(output.window1.maxidx(1),output.window1.maxidx(2)) = NaN;
% dataMatWindow(output.window1.minidx(1),output.window1.minidx(2)) = NaN;
% contourf(dataMatWindow,'edgecolor','none')

%% determine window2 range around window1 min/max activity

% default window2 range is the range/size of window1, so window1 and window2 should give same results if window2 range is not specified

if ~isfield(cfg,'window2_xrange')
    cfg.window2_xrange = cfg.window1_xlim(end) - cfg.window1_xlim(1);
end

if ~isfield(cfg,'window2_yrange') 
    cfg.window2_yrange = cfg.window1_ylim(end) - cfg.window1_ylim(1);
end

%% window2 around maximum window1 activity    

output.window2.measure1 = 'activity around max';
% frequency range
tempY = [ceil(output.window1.maxfoi - cfg.window2_yrange/2) floor(output.window1.maxfoi + cfg.window2_yrange/2)]; % determine y min and max values (freq min/max)

% ensure window2 is inside window1
if strcmpi(cfg.win2insidewin1,'yes')
    if tempY(1) < cfg.window1_ylim(1)
        tempY = cfg.window1_ylim(1) - tempY(1) + tempY;
    end
    if tempY(2) > cfg.window1_ylim(2)
        tempY = cfg.window1_ylim(2) - tempY(2) + tempY;
    end
end

output.window2.yrange_maxidx = dsearchn(cfg.yfoi',tempY'); % find matching indices
output.window2.yrange_maxfoi = [cfg.yfoi(output.window2.yrange_maxidx(1)) cfg.yfoi(output.window2.yrange_maxidx(2))]; % find matching frequencies

% time range
tempX = [output.window1.maxtoi - cfg.window2_xrange/2 output.window1.maxtoi + cfg.window2_xrange/2];

% ensure window2 is inside window1
if strcmpi(cfg.win2insidewin1,'yes')
    if tempX(1) < cfg.window1_xlim(1)
        tempX = cfg.window1_xlim(1) - tempX(1) + tempX;
    end
    if tempX(2) > cfg.window1_xlim(2)
        tempX = cfg.window1_xlim(2) - tempX(2) + tempX;
    end
end

output.window2.xrange_maxidx = dsearchn(cfg.xtoi',tempX');
output.window2.xrange_maxtoi = [cfg.xtoi(output.window2.xrange_maxidx(1)) cfg.xtoi(output.window2.xrange_maxidx(2))];

% select data for window2
window2data = dataMat;
% x: convert values outside window to NaN
tempxlimidx = cfg.xtoi >= output.window2.xrange_maxtoi(1) & cfg.xtoi <= output.window2.xrange_maxtoi(2); % indices inside window
window2data(:,~tempxlimidx) = NaN;
% y: convert values outside window to NaN
tempylimidx = cfg.yfoi >= output.window2.yrange_maxfoi(1) & cfg.yfoi <= output.window2.yrange_maxfoi(2); % indices inside window
window2data(~tempylimidx,:) = NaN;
% compute mean in window2
output.window2.max_mean = nanmean(window2data(:));

% create thresholded map for plotting
window2data(~isnan(window2data)) = 1;
window2data(isnan(window2data)) = 0;
output.window2.max_thresholdmap = window2data;
% rectangle plotting position
output.window2.rectangleposition_max = [output.window2.xrange_maxtoi(1) output.window2.yrange_maxfoi(1) output.window2.xrange_maxtoi(2)-output.window2.xrange_maxtoi(1) output.window2.yrange_maxfoi(2)-output.window2.yrange_maxfoi(1)];
    
% figure(2), clf
% contourf(window2data,'edgecolor','none')

%% window2 around minimum window1 activity    

output.window2.measure2 = 'activity around min';
% frequency range
tempY = [ceil(output.window1.minfoi - cfg.window2_yrange/2) floor(output.window1.minfoi + cfg.window2_yrange/2)]; % determine y min and max values (freq min/max)

% ensure window2 is inside window1
if strcmpi(cfg.win2insidewin1,'yes')
    if tempY(1) < cfg.window1_ylim(1)
        tempY = cfg.window1_ylim(1) - tempY(1) + tempY;
    end
    if tempY(2) > cfg.window1_ylim(2)
        tempY = cfg.window1_ylim(2) - tempY(2) + tempY;
    end
end

output.window2.yrange_minidx = dsearchn(cfg.yfoi',tempY'); % find matching indices
output.window2.yrange_minfoi = [cfg.yfoi(output.window2.yrange_minidx(1)) cfg.yfoi(output.window2.yrange_minidx(2))]; % find matching frequencies

% time range
tempX = [output.window1.mintoi - cfg.window2_xrange/2 output.window1.mintoi + cfg.window2_xrange/2];

% ensure window2 is inside window1
if strcmpi(cfg.win2insidewin1,'yes')
    if tempX(1) < cfg.window1_xlim(1)
        tempX = cfg.window1_xlim(1) - tempX(1) + tempX;
    end
    if tempX(2) > cfg.window1_xlim(2)
        tempX = cfg.window1_xlim(2) - tempX(2) + tempX;
    end
end

output.window2.xrange_minidx = dsearchn(cfg.xtoi',tempX');
output.window2.xrange_mintoi = [cfg.xtoi(output.window2.xrange_minidx(1)) cfg.xtoi(output.window2.xrange_minidx(2))];

% select data for window2
window2data = dataMat;
% x: convert values outside window to NaN
tempxlimidx = cfg.xtoi >= output.window2.xrange_mintoi(1) & cfg.xtoi <= output.window2.xrange_mintoi(2); % indices inside window
window2data(:,~tempxlimidx) = NaN;
% y: convert values outside window to NaN
tempylimidx = cfg.yfoi >= output.window2.yrange_minfoi(1) & cfg.yfoi <= output.window2.yrange_minfoi(2); % indices inside window
window2data(~tempylimidx,:) = NaN;
% compute mean in window2
output.window2.min_mean = nanmean(window2data(:));

% create thresholded map for plotting
window2data(~isnan(window2data)) = 1;
window2data(isnan(window2data)) = 0;
output.window2.min_thresholdmap = window2data;
% rectangle plotting position
output.window2.rectangleposition_min = [output.window2.xrange_mintoi(1) output.window2.yrange_minfoi(1) output.window2.xrange_mintoi(2)-output.window2.xrange_mintoi(1) output.window2.yrange_minfoi(2)-output.window2.yrange_minfoi(1)];
    
% figure(2), clf
% contourf(window2data,'edgecolor','none')

%% plot

if ~isfield(cfg,'plotstat') cfg.plotstat = 'no'; end
if ismember(cfg.plotstat,{'min','max','minmax','maxmin','yes'})
    tempcfg = [];
    tempcfg.xtoi = cfg.xtoi;
    tempcfg.yfoi = cfg.yfoi;
    
    if ~isfield(cfg,'title')
        tempcfg.title = cfg.chan;
    else
        tempcfg.title = {cfg.title cfg.chan};
    end 
    
    if ~isfield(cfg,'fontsize')
        tempcfg.fontsize = 9;
    else
        tempcfg.fontsize = cfg.fontsize;
    end
    
    if ~isfield(cfg,'clim')
        tempcfg.clim = [-2.5 2.5];
    else
        tempcfg.clim = cfg.clim;
    end

    if ~isfield(cfg,'ylog') cfg.ylog = 'yes'; end
    tempcfg.ylog = cfg.ylog;
    if isfield(cfg,'subplot')
        tempcfg.subplot = cfg.subplot;
    end
    if isfield(cfg,'saveas')
        tempcfg.saveas = cfg.saveas;
    end
    
    tempout = ft_contourfTFR(tempcfg,dataMat); % plot data
    output.plotcfg = tempout;
    
    % highlight max/min values
    hold on
    plot(output.window1.maxtoi,output.window1.maxfoi,'r*','markersize',7)
    plot(output.window1.mintoi,output.window1.minfoi,'m*','markersize',7)
   
    % draw rectangles
    rectangle('Position',output.window1.rectangleposition,'edgecolor','w','linewidth',2,'linestyle',':') % draw input window (selection window)     
    if strcmpi(cfg.plotstat,'min')
        rectangle('Position',output.window2.rectangleposition_min,'edgecolor','m','linewidth',2,'linestyle','-') % draw output window (window centered around min value)
    elseif strcmpi(cfg.plotstat,'max')
        rectangle('Position',output.window2.rectangleposition_max,'edgecolor','r','linewidth',2,'linestyle','-') % draw output window (window centered around max value)
    else
        rectangle('Position',output.window2.rectangleposition_min,'edgecolor','m','linewidth',2,'linestyle','-') % draw output window (window centered around min value)
        rectangle('Position',output.window2.rectangleposition_max,'edgecolor','r','linewidth',2,'linestyle','-') % draw output window (window centered around max value)
    end
    
end

output.cfg = cfg;

end