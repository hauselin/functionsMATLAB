function [dataout] = ft_computeWITPC(cfg,data)

%% check input data and organize

if ~isfield(cfg,'designmatrix') 
    error('Provide design matrix!'); 
end

if ~isfield(data,'fourierspctrm') % see if fourier coefficients exist
    error('Make sure data contain fourierspctrm field!');
end

if ~isfield(cfg,'trials')
    cfg.trials = 'all';
end

if size(cfg.designmatrix,2) > size(cfg.designmatrix,1) % tranpose design matrix such that it's trial_regressor
    cfg.designmatrix = cfg.designmatrix';
end

if size(cfg.designmatrix,1) ~= size(data.fourierspctrm,1) % check if number of trials in design matrix and data match
    if isfield(cfg,'trials') % if selecting trials
        data.fourierspctrm = data.fourierspctrm(cfg.trials,:,:,:);
        data.cumtapcnt = data.cumtapcnt(cfg.trials,:);
    else
        error('Trial numbers in design matrix and data do not match!');
    end
end

if ~isfield(cfg,'n_permutes') cfg.n_permutes = 0; end

%% select channels and trials

if isfield(cfg,'chan')
    tempcfg = [];
    tempcfg.channel = cfg.chan;
    tempcfg.trials = cfg.trials;
    tempcfg.avgoverchan = 'no';
    tempcfg.avgoverrpt = 'no';
    data = ft_selectdata(tempcfg,data);
end

%% prepare output structure

tempcfg = [];
tempcfg.avgoverrpt = 'yes';
dataout = data;
dataout.fourierspctrm = dataout.fourierspctrm(1,:,:,:);
dataout = ft_selectdata(tempcfg,data);
    
%% calculate phase angle

data.fourierspctrm = angle(data.fourierspctrm); % calculate phase angle

%% calculate ITPC 

dataout.itpc = abs( mean(exp(1i*squeeze(data.fourierspctrm)),1) ); % output: chan_freq_time

%% calculate weighted ITPC: phase modulation by trial-varying information

dataout.witpc = dataout.fourierspctrm;
if cfg.n_permutes > 0
    dataout.witpc_z = dataout.fourierspctrm; % permuted witpc zscores
else
    dataout.witpc_z = [];
end

ntrials = length(cfg.designmatrix);
modulator = cfg.designmatrix;

% chanI = 2
for chanI = 1:length(dataout.label)

    % data for this channel
    tempData = squeeze(data.fourierspctrm(:,chanI,:,:));

    % observed weighted ITPC for this channel
    realwitpc = squeeze( abs( mean(modulator .* exp(1i*tempData),1) ) ); % rts modulating length of phase angles
    dataout.witpc(chanI,:,:) = realwitpc;

    % permutation testing
    if cfg.n_permutes > 0
        perm_witpc = zeros(cfg.n_permutes,size(realwitpc,1),size(realwitpc,2));
        for permI = 1:cfg.n_permutes
            perm_witpc(permI,:,:) = squeeze( abs( mean(modulator(randperm(ntrials)) .* exp(1i*tempData),1) ) ); % rts modulating length of phase angles
            disp(['Generated permuted distribution ' num2str(permI) ' of ' num2str(cfg.n_permutes) ' (channel ' num2str(chanI) ' of ' num2str(length(data.label)) ')']);
        end

        witpc_z = (realwitpc - squeeze(nanmean(perm_witpc,1))) ./ squeeze(nanstd(perm_witpc,[],1));
        dataout.witpc_z(chanI,:,:) = witpc_z;
    end
    
end

%% prepare output

dataout = rmfield(dataout,'fourierspctrm');
dataout.cfg = cfg;

end