function [dataout] = ft_regressTFRPower_permute(cfg,data)
% Performs multiple regression, permutation, and cluster correction.
%
% Written in MATLAB R2017b
% Last modified by Hause Lin 21-08-18 20:22 hauselin@gmail.com

%% Check input 

if ~isfield(cfg,'designmatrix') 
    error('Provide design matrix!'); 
end

if size(cfg.designmatrix,1) ~= size(data.powspctrm,1) % check if number of trials in design matrix and data match
    if isfield(cfg,'trials') % if selecting trials
        data.powspctrm = data.powspctrm(cfg.trials,:,:,:);
        data.cumtapcnt = data.cumtapcnt(cfg.trials,:);
    else
        error('Trial numbers in design matrix and data do not match!');
    end
end

%% Set default values

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
else
    designMatrix = cfg.designmatrix; % don't add intercept/constant
end

%% Create structure to store output

dataout = data;
dataout.powspctrm = NaN(size(designMatrix,2),size(data.powspctrm,2),size(data.powspctrm,3),size(data.powspctrm,4)); % permuted zmap (regressorZ_chan_freq_time)
dataout.cumtapcnt = data.cumtapcnt(1:size(designMatrix,2),:); 

if strcmpi(cfg.type,'zpermutecluster') % save max cluster information (cluster correction)
    dataout.cfg.maxclusterinfo = NaN(size(data.powspctrm,2),cfg.n_permutes,size(designMatrix,2)); 
end

%% Compute real beta values, zmap, and cluster correction

for chanI = 1:length(data.label) % for each channel, perform regression and permutation to get zmap
    
    % reshape data for single-trial regression
    tf3d = squeeze(data.powspctrm(:,chanI,:,:)); % trial_freq_time
    tf2d = reshape(tf3d,size(designMatrix,1),length(data.time)*length(data.freq)); % reshape to 2d (trial_freq*time)
    
    % ranktransform power if necesary
    if strcmpi(cfg.transformY,'tiedrank')
        tf2d = tiedrank(tf2d);
    end
    
    % regression using least-squares equation (real betas)
    betaCoefs = (designMatrix'*designMatrix)\designMatrix'*tf2d; % designMat\tf2d; regression (least-squares equation))
    % size(betaCoefs)
    
    % reshape real beta matrix
    realbeta = reshape(betaCoefs,[size(designMatrix,2) length(data.freq) length(data.time)]); % reshape to 3d matrix 
    
    disp(['Computed real beta values with ordinary least-squares regression (channel ' num2str(chanI) ' of ' num2str(length(data.label)) ')']);
    
    if strcmpi(cfg.type,'ols') % if ols, save realbeta values and move to next loop/channel
        dataout.powspctrm(:,chanI,:,:) = realbeta;
        continue
    end
    
    %% permute to get permuted distributions
    
    % initialize null hypothesis matrices
    permuted_bvals = zeros(cfg.n_permutes,size(designMatrix,2),length(data.freq),length(data.time)); % permutes_betas_frex_time
    % size(permuted_bvals)
    max_clust_info = zeros(cfg.n_permutes,size(designMatrix,2)); % permute_betas
    % size(max_clust_info)
    
    % generate pixel-specific null hypothesis parameter distributions
    for permi = 1:cfg.n_permutes
        % randomly shuffle trial order
        fakedesignMat = designMatrix(randperm(size(designMatrix,1)),:);

        % compute beta-map of null hypothesis
        fakebeta = (fakedesignMat'*fakedesignMat)\fakedesignMat'*tf2d;

        % reshape to 2D map for cluster-correction
        fakebeta = reshape(fakebeta,[size(designMatrix,2) length(data.freq) length(data.time)]); % reshape to 3d matrix (coefficient_freq_time)
        
        % save all permuted values
        permuted_bvals(permi,:,:,:) = fakebeta;
        
        disp(['Generated permuted distribution ' num2str(permi) ' of ' num2str(cfg.n_permutes) ' (channel ' num2str(chanI) ' of ' num2str(length(data.label)) ')']);
    end

    %% compute zmap: (realbeta / mean of permuted values) / sd of permuted values
    
    zmap = NaN([size(designMatrix,2) length(data.freq) length(data.time)]);
    for betaI = 1:size(designMatrix,2) % for each beta, compute permuted z value
        zmap(betaI,:,:) = (squeeze(realbeta(betaI,:,:))-squeeze(nanmean(permuted_bvals(:,betaI,:,:),1))) ./ squeeze(nanstd(permuted_bvals(:,betaI,:,:),[],1));
        disp(['Computed permuted zmap regressor ' num2str(betaI) ' of ' num2str(size(designMatrix,2)) ' (channel ' num2str(chanI) ' of ' num2str(length(data.label)) ')']);
    end
    
    dataout.powspctrm(:,chanI,:,:) = zmap;
    
    if strcmpi(cfg.type,'zpermute') % if zpermute, save zmap values and move to next loop/channel
        dataout.powspctrm(:,chanI,:,:) = zmap;
        continue
    end
    
    %% cluster correction on permuted data
    
    for permi = 1:cfg.n_permutes % for each permuted map
        for betaI = 1:size(designMatrix,2) % for each beta coefficient
            % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
            fakecorrsz = squeeze( (permuted_bvals(permi,betaI,:,:)-nanmean(permuted_bvals(:,betaI,:,:),1) ) ./ nanstd(permuted_bvals(:,betaI,:,:),[],1) );
            fakecorrsz(abs(fakecorrsz) < norminv(1-cfg.voxel_pvalue)) = 0; % identify values greater/smaller than voxel_pvalue
            clustinfo = bwconncomp(fakecorrsz); % get number of elements in largest supra-threshold cluster
            max_clust_info(permi,betaI) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
            
        end
        disp(['Cluster correction ' num2str(permi) ' of ' num2str(cfg.n_permutes) ' (channel ' num2str(chanI) ' of ' num2str(length(data.label)) ')']);
    end
    
    dataout.cfg.maxclusterinfo(chanI,:,:) = max_clust_info;
    disp(['Finished cluster correction (channel ' num2str(chanI) ' of ' num2str(length(data.label)) ')']);
    
end

%% format output data structure

if size(dataout.powspctrm,1) == 1 % if there's only one coefficient, squeeze it
    dataout.powspctrm = squeeze(dataout.powspctrm);
end

dataout.cfg.regressparams = cfg; % save regression parameters in ft structure

end