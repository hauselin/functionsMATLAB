function [dataout] = ft_regressERF_permute(cfg,data)
% Performs multiple regression, permutation, and cluster correction.
%
% Written in MATLAB R2017b
% Last modified by Hause Lin 22-08-18 13:09 hauselin@gmail.com

%% Check input 

if ~isfield(cfg,'designmatrix') 
    error('Provide design matrix!'); 
end

if isfield(data,'cfg') % if data is in fieldtrip format, convert back to EEG.data format
    % do shit here
end

if size(cfg.designmatrix,1) ~= size(data.data,3) % check if number of trials in design matrix and data match
    if isfield(cfg,'trials')
        data = pop_select(data,'trial',cfg.trials); % select only clean epochs 
        if size(cfg.designmatrix,1) ~= size(data.data,3)
            error('Trial numbers in design matrix and data do not match!');
        end
    else
        error('Trial numbers in design matrix and data do not match!');
    end
end

%% Set default values

if ~isfield(cfg,'chan') % default all channels
    cfg.chan = {data.chanlocs.labels};
    Y = data.data;
    cfg.chanind = 1:length({data.chanlocs.labels});
else % select channels
    if isa(cfg.chan,'char')
        cfg.chan = {cfg.chan};
    end
    cfg.chanind = [];
    for chani = 1:length(cfg.chan)
       cfg.chanind(chani) = find(strcmpi({data.chanlocs.labels},cfg.chan{chani}));
    end
    Y = data.data(cfg.chanind,:,:);
end

if ~isfield(cfg,'times') % default all time points (in seconds)
    cfg.times = data.times / 1000;
    cfg.timesind = 1:length(data.times);
else
    times = dsearchn(data.times'/1000,[cfg.times(1) cfg.times(end)]');
    cfg.times = data.times(times(1):times(2))/1000;
    cfg.timesind = times(1):times(end);
    Y = squeeze(Y(:,cfg.timesind,:));
end

cfg.trials = data.trials;

cfg.datadim = [size(Y,1) size(Y,2) size(Y,3)];

if ~isfield(cfg,'intercept') cfg.intercept = 'no'; end
if ~isfield(cfg,'model') cfg.model = 'regression'; end
if ~isfield(cfg,'transformY') cfg.transformY = 'no'; end
if ~isfield(cfg,'type') cfg.type = 'ols'; end % ols, zpermute, zpermutecluster
if ~isfield(cfg,'voxel_pvalue') cfg.voxel_pvalue = 0.05; end
if ~isfield(cfg,'mcc_cluster_pvalue') cfg.mcc_cluster_pvalue = 0.05; end
if ~isfield(cfg,'n_permutes') cfg.n_permutes = 500; end

if size(cfg.designmatrix,2) > size(cfg.designmatrix,1) % tranpose design matrix such that it's trial_regressor
    cfg.designmatrix = cfg.designmatrix';
end

%% Add intercept if necessary

if strcmpi(cfg.intercept,'yes')
    designMatrix = [ones(size(cfg.designmatrix,1),1) cfg.designmatrix]; % add intercept/constant
    regressors = [{'b0'} cfg.regressors];
else
    designMatrix = cfg.designmatrix; % don't add intercept/constant
    regressors = cfg.regressors;
end

%% Create structure to store output

dataout = struct();

%% Compute real beta values, zmap, and cluster correction

% get data and reshape data to trial_datapoints
data2d = reshape(Y,length(cfg.chan)*length(cfg.times),cfg.trials)'; % trial_datapoints

% ranktransform data if necesary
if strcmpi(cfg.transformY,'tiedrank')
    data2d = tiedrank(data2d);
end

% regression using least-squares equation (real betas)
betaCoefs = (designMatrix'*designMatrix)\designMatrix'*data2d; % designMat\data2d; regression (least-squares equation))
% size(betaCoefs)

% reshape real beta matrix
realbeta = reshape(betaCoefs,[size(designMatrix,2) length(cfg.chan) length(cfg.times)]); % reshape to 3d matrix 

disp('Computed real beta values with ordinary least-squares regression');

dataout.realbeta = realbeta;

%% permute to get permuted distributions
    
if strcmpi(cfg.type,'zpermute') % if zpermute
    
    % initialize null hypothesis matrices
    permuted_bvals = zeros(cfg.n_permutes,size(designMatrix,2),length(cfg.chan),length(cfg.times)); % permutes_betas_chan_time
    % size(permuted_bvals)
    
    % generate pixel-specific null hypothesis parameter distributions
    for permi = 1:cfg.n_permutes
        % randomly shuffle trial order
        fakedesignMat = designMatrix(randperm(size(designMatrix,1)),:);

        % compute beta-map of null hypothesis
        fakebeta = (fakedesignMat'*fakedesignMat)\fakedesignMat'*data2d;

        % reshape to 2D map for cluster-correction
        fakebeta = reshape(fakebeta,[size(designMatrix,2) length(cfg.chan) length(cfg.times)]); % reshape to 3d matrix, coefficient_chan_time

        % save all permuted values
        permuted_bvals(permi,:,:,:) = fakebeta;
        
        disp(['Generated permuted distribution ' num2str(permi) ' of ' num2str(cfg.n_permutes)]);
    end

    %% compute zmap: (realbeta / mean of permuted values) / sd of permuted values
    
    zmap = NaN([size(designMatrix,2) length(cfg.chan) length(cfg.times)]);
    for betaI = 1:size(designMatrix,2) % for each beta, compute permuted z value
        zmap(betaI,:,:) = (squeeze(realbeta(betaI,:,:))-squeeze(nanmean(permuted_bvals(:,betaI,:,:),1))) ./ squeeze(nanstd(permuted_bvals(:,betaI,:,:),[],1));
        disp(['Computed permuted zmap regressor ' num2str(betaI) ' of ' num2str(size(designMatrix,2))]);
    end
    
    dataout.zmap = zmap;
end

%% cluster correction on permuted data

if strcmpi(cfg.type,'zpermutecluster') % if zpermutecluster
    
    max_clust_info = zeros(cfg.n_permutes,size(designMatrix,2)); % permute_betas
    % size(max_clust_info)
    
    for permi = 1:cfg.n_permutes % for each permuted map
        for betaI = 1:size(designMatrix,2) % for each beta coefficient
            % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
            fakecorrsz = squeeze( (permuted_bvals(permi,betaI,:,:)-nanmean(permuted_bvals(:,betaI,:,:),1) ) ./ nanstd(permuted_bvals(:,betaI,:,:),[],1) );
            fakecorrsz(abs(fakecorrsz) < norminv(1-cfg.voxel_pvalue)) = 0; % identify values greater/smaller than voxel_pvalue
            clustinfo = bwconncomp(fakecorrsz); % get number of elements in largest supra-threshold cluster
            max_clust_info(permi,betaI) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
        end
        disp(['Cluster correction ' num2str(permi) ' of ' num2str(cfg.n_permutes)]);
    end
    
    dataout.maxclusterinfo = max_clust_info;
end

%% organise output structure (suitable for fieldtrip processing)

dataft = eeglab2fieldtrip(data,'preprocessing');

elec = dataft.elec;
elec.pnt = elec.pnt(cfg.chanind,:);
elec.label = elec.label(1,cfg.chanind);

dataout.cfg = cfg; % save regression parameters in structure
dataout.regressors = regressors;
dataout.designmat = designMatrix;
dataout.time = cfg.times;
dataout.dimord = 'chan_time';
dataout.label = cfg.chan';
dataout.elec = elec;
dataout.chanloc = data.chanlocs(cfg.chanind);
dataout.cfg.subject = data.subject;

end