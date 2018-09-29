function [cfg] = ft_contourfTFR(cfg,data)
% Plots time-frequency map.
% Data should be 2D (time_freq matrix) or the output from ft_freqanalysis.
% If data is the output from ft_freqanalysis, the function will
% try to plot data from FCz channel; but if FCz is unavailable,
% it will plot the first channel.
%
% cfg.xtoi: time in seconds (not necessary if input is ft structure)
% cfg.yfoi: frequencies (not necessary if input is ft structure)
% cfg.ylog = 'yes' (default) or 'no'
% cfg.subplot: subplot number (default [1 1 1]);
% cfg.parameter: field/parameter to plot (default 'powspctrm');
% cfg.chan: which channel to plot (if input is ft structure: default is FCz; if not, just plot channel 1); 
% cfg.trials: which trial of the 'rpt' dimension to plot (if input is ft structure; default trial 1)
% cfg.xlim: x (time) limits for plot (default [-0.9 0.9])
% cfg.ylim: y (freq) limits for plot (default [1 max(cfg.yfoi)])
% cfg.yticksn: number of y-axis tick marks (default 13)
% cfg.ylog: log10 y-axis (default 'yes')
% cfg.colormap: 'viridis' (default; default is also 'parula')
% cfg.clim: colour limits (default [-2 2])
% cfg.square = 'yes' (default) or 'no'
% cfg.colorbar = 'yes' (default) or 'no'
% cfg.colorbarticks: 'yes' (default) or 'no'
% cfg.xlabel: 'Time (s)'
% cfg.ylabel: 'Frequency (Hz)'
% cfg.title: plot title (default: [cfg.chan ' power' ])
% cfg.fontsize: 14 (default)
% cfg.saveas: specify where to save and file name
% 
% Written in MATLAB R2017b
% Last modified by Hause Lin 05-09-18 21:50 hauselin@gmail.com

%% fix structure if input data isn't powspctrm

if ~isfield(cfg,'parameter')
    cfg.parameter = 'powspctrm';
else
    data.powspctrm = data.(cfg.parameter);
end

%% check input data structure

if ~isstruct(data) && ndims(data) == 2 && ismatrix(data) % if 2D matrix, then just plot the data (need to provide time and freq: xtoi and yfoi)
    if ~isfield(cfg,'chan') cfg.chan = 'Channel'; end
    if ~isfield(cfg,'xtoi') error('Specify cfg.xtoi'); end
    if ~isfield(cfg,'yfoi') error('Specify cfg.yfoi'); end
    dataMat = data;
elseif isstruct(data) && ndims(data.powspctrm) == 3 % multi-channel data
    if ~isfield(cfg,'chan') cfg.chan = 'FCz'; end
    cfg.chanIndex = find(strcmpi(cfg.chan,data.label));
    if isempty(cfg.chanIndex) 
        cfg.chanIndex = 1; 
        cfg.chan = data.label{cfg.chanIndex};
    end
    cfg.xtoi = data.time;
    cfg.yfoi = data.freq;
    dataMat = squeeze(data.powspctrm(cfg.chanIndex,:,:));
elseif isstruct(data) && ndims(data.powspctrm) == 4 % multi-trial and channel data
    if ~isfield(cfg,'chan') cfg.chan = 'FCz'; end
    cfg.chanIndex = find(strcmpi(cfg.chan,data.label));
    if isempty(cfg.chanIndex) 
        cfg.chanIndex = 1; 
        cfg.chan = data.label{cfg.chanIndex};
    end
    cfg.xtoi = data.time;
    cfg.yfoi = data.freq;
    if ~isfield(cfg,'trials') cfg.trials = 1:size(data.powspctrm,1); end % default plot mean of all trials
    dataMat = squeeze(nanmean(data.powspctrm(cfg.trials,cfg.chanIndex,:,:),1));
else
    error('Check input data!')
end

if ~isfield(cfg,'xlim') cfg.xlim = [-0.9 0.9]; end
if ~isfield(cfg,'ylim') cfg.ylim = [1 max(cfg.yfoi)]; end
if ~isfield(cfg,'yticksn') cfg.yticksn = 13; end
if ~isfield(cfg,'ylog') cfg.ylog = 'yes'; end
if strcmpi(cfg.ylog,'yes')
    cfg.yticks = round(logspace(log10(cfg.yfoi(1)),log10(cfg.yfoi(end)),cfg.yticksn),1);
else
    cfg.yticks = round(linspace(cfg.yfoi(1),cfg.yfoi(end),cfg.yticksn));
end
if ~isfield(cfg,'colormap') cfg.colormap = 'viridis'; end
if ~isfield(cfg,'clim') 
    if strcmpi(cfg.parameter,'itpc')
        cfg.clim = [0 0.4];
    else
        cfg.clim = [-2 2];
    end
end
if ~isfield(cfg,'square') cfg.square = 'yes'; end
if ~isfield(cfg,'colorbar') cfg.colorbar = 'yes'; end
if ~isfield(cfg,'colorbarticks') cfg.colorbarticks = 'yes'; end
if ~isfield(cfg,'xlabel') cfg.xlabel = 'Time (s)'; end
if ~isfield(cfg,'ylabel') cfg.ylabel = 'Frequency (Hz)'; end
if ~isfield(cfg,'title') cfg.title = cfg.chan; end
if ~isfield(cfg,'fontsize') cfg.fontsize = 14; end

% plot
if ~isfield(cfg,'subplot') 
    figure
else
    subplot(cfg.subplot(1),cfg.subplot(2),cfg.subplot(3))
end

% set color map
try 
    colormap(cfg.colormap)
catch
    try 
        colormap('viridis')
    catch
        try
            colormap('parula')
        catch
            colormap('jet')
        end
    end
end

%% actual plotting 

contourf(cfg.xtoi,cfg.yfoi,dataMat,40,'linecolor','none')
if strcmpi(cfg.ylog,'yes')
    set(gca,'xlim',cfg.xlim,'clim',cfg.clim,'yscale','log','ytick',cfg.yticks,'FontSize',cfg.fontsize,'yminortick','off')
else
    set(gca,'xlim',cfg.xlim,'clim',cfg.clim,'ytick',cfg.yticks,'FontSize',cfg.fontsize)
end
xlabel(cfg.xlabel); ylabel(cfg.ylabel); title(cfg.title);
if strcmpi(cfg.colorbar,'yes') cb = colorbar; end

if isfield(cfg,'square')
    axis square
end

if isfield(cfg,'saveas')
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
    print(gcf,'-djpeg','-r200',cfg.saveas);
end

end