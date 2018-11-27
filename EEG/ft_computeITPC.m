function [dataout] = ft_computeITPC(cfg,data)
% Computes ITPC values from single-trial fourier coefficients (output from ft_freqanalysis) 
% 
% cfg.trials: which trials to use (default 'all')
% cfg.avgoverrpt: average trials or not (default 'yes')
% cfg.itpcz: transform ITPC to ITPCz (deafult 'no')
% 
% Usage
% cfg = [];
% ft_computeITPC(cfg,data)
%
% Written by Hause Lin 14-11-18 21:31 hauselin@gmail.com

%% check input data

if ~isfield(data,'fourierspctrm') % see if fourier coefficients exist
    error('Check input data');
else
    data.fourierspctrm = angle(data.fourierspctrm); % get phase angle
end

if ~isfield(cfg,'trials') % select trials
    cfg.trials = 'all';
end

if ~isfield(cfg,'avgoverrpt') % average over trials
    cfg.avgoverrpt = 'yes';
end

if ~isfield(cfg,'itpcz') % transform ITPC to ITPCz
    cfg.itpcz = 'no';
end

%% prepare output structure

if strcmpi(cfg.avgoverrpt,'yes') % return average
    tempcfg = [];
    tempcfg.avgoverrpt = cfg.avgoverrpt;
    dataout = data;
    dataout.fourierspctrm = dataout.fourierspctrm(1,:,:,:);
    dataout = ft_selectdata(tempcfg,data);
else % return single-trial data
    tempcfg = [];
    tempcfg.trials = cfg.trials;
    tempcfg.avgoverrpt = cfg.avgoverrpt;
    tempcfg.avgoverrpt = 'no';
    dataout = ft_selectdata(tempcfg,data);
end
    
%% calculate ITPC 

% select trials
if ~strcmpi(cfg.trials,'all') % itpc on select trials
    dataSubset = data.fourierspctrm(cfg.trials,:,:,:);
    disp(['Averaging ' num2str(length(cfg.trials)) ' trials']);
else
    dataSubset = data.fourierspctrm;
end

if strcmpi(cfg.avgoverrpt,'yes') % return average
    dataout.fourierspctrm = squeeze( abs( mean(exp(1i*dataSubset),1) ) ); % output: chan_freq_time
else % return single-trial data
    dataout.fourierspctrm = exp(1i*dataSubset); % output: trial_chan_freq_time
    cfg.info = 'abs(mean(x.itpc,1)) to get mean itpc from single trial data';
end

if strcmpi(cfg.itpcz,'yes') % transform ITPC to ITPCz (see Cohen, 2014, p.249)
    dataout.fourierspctrm = length(cfg.trials) .* (dataout.fourierspctrm.^2);
end

%% prepare output

dataout = renamefields(dataout,'fourierspctrm','itpc');
dataout.cfg = cfg;

end