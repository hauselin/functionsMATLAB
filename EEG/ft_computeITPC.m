function [dataout] = ft_computeITPC(cfg,data)

%% check input data

if ~isfield(data,'fourierspctrm') % see if fourier coefficients exist
    error('Check input data');
else
    data.fourierspctrm = angle(data.fourierspctrm); % get phase angle
end

if ~isfield(cfg,'trials')
    cfg.trials = 'all';
end

if ~isfield(cfg,'avgoverrpt')
    cfg.avgoverrpt = 'yes';
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
else
    dataSubset = data.fourierspctrm;
end

if strcmpi(cfg.avgoverrpt,'yes') % return average
    dataout.fourierspctrm = squeeze( abs( mean(exp(1i*dataSubset),1) ) ); % output: chan_freq_time
else % return single-trial data
    dataout.fourierspctrm = exp(1i*dataSubset); % output: trial_chan_freq_time
    cfg.info = 'abs(mean(x.itpc,1)) to get mean itpc from single trial data';
end

if isfield(cfg,'itpc_rayleighZ') % transform ITPC to ITPCz (see Cohen, 2014, p.249)
    dataout.fourierspctrm = length(cfg.trials) .* (dataout.fourierspctrm.^2);
end

%% prepare output

dataout = renamefields(dataout,'fourierspctrm','itpc');
dataout.cfg = cfg;

end