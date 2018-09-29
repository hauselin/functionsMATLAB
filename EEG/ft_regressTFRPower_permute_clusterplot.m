function [dataout] = ft_regressTFRPower_permute_clusterplot(cfg,data)
%
% Written in MATLAB R2017b
% Last modified by Hause Lin 21-08-18 20:22 hauselin@gmail.com

%% Check input 

if ~isfield(cfg,'chan') 
    error('Indicate one channel!'); 
end

if isa(cfg.chan,'cell')
    error('Indicate only one channel!'); 
end

if ~isfield(data.cfg,'maxclusterinfo') 
    error('No cluster correction information! Rerun ft_regressTFRPower_permute.'); 
end

%% Set default values

if ~isfield(cfg,'regressors') cfg.regressors = 1:size(data.powspctrm,1); end
if ~isfield(cfg,'voxel_pvalue') cfg.voxel_pvalue = 0.05; end
if ~isfield(cfg,'mcc_cluster_pvalue') cfg.mcc_cluster_pvalue = 0.05; end
if ~isfield(cfg,'clim') cfg.clim = [-2 2]; end

%% Get channel data

zmap = squeeze(data.powspctrm(:,find(ismember(data.label,cfg.chan)),:,:));
max_clust_info = squeeze(data.cfg.maxclusterinfo(find(ismember(data.label,cfg.chan)),:,:));

%% Plot

fig = figure(183249479);
set(fig,'NumberTitle','off','Name','z score cluster correction'); clf
set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure

for regressorI=1:length(cfg.regressors)
    
    colormap viridis
    
    %% get data
    
    zmapI = squeeze(zmap(cfg.regressors(regressorI),:,:));
    
    %% permuted z map
    
    subplot(3,length(cfg.regressors),regressorI)
    contourf(data.cfg.toi,data.cfg.foi,zmapI,40,'linecolor','none')
    axis square
    set(gca,'clim',[cfg.clim(1) cfg.clim(2)])
    if isfield(cfg,'regressorlabels')
        title({cfg.regressorlabels{regressorI} 'Unthresholded z map'})
    end
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
     
    %% apply uncorrected threshold
    
    zmapthresh = zmapI;
    zmapthresh(abs(zmapthresh) < norminv(1-cfg.voxel_pvalue)) = 0;
    subplot(3,length(cfg.regressors),regressorI+length(cfg.regressors))
    contourf(data.cfg.toi,data.cfg.foi,zmapthresh,40,'linecolor','none')
    axis square
    set(gca,'clim',[cfg.clim(1) cfg.clim(2)])
    title('Uncorrected thresholded z map')
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    
    %% apply cluster-level corrected threshold
    
    zmapthresh = zmapI;
    % uncorrected pixel-level threshold
    zmapthresh(abs(zmapthresh) < norminv(1-cfg.voxel_pvalue)) = 0;
    % find islands and remove those smaller than cluster size threshold
    clustinfo = bwconncomp(zmapthresh);
    clust_info = cellfun(@numel,clustinfo.PixelIdxList); % number of elements in each cluster
    clust_threshold = prctile(max_clust_info(:,regressorI),100-cfg.mcc_cluster_pvalue*100);
    
    % identify clusters to remove
    whichclusters2remove = find(clust_info < clust_threshold);
    
    % remove clusters
    for i=1:length(whichclusters2remove)
        zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)}) = 0;
    end
    
    subplot(3,length(cfg.regressors),regressorI+length(cfg.regressors)*2)
    contourf(data.cfg.toi,data.cfg.foi,zmapthresh,40,'linecolor','none')
    axis square
    set(gca,'clim',[cfg.clim(1) cfg.clim(2)])
    title('Cluster-corrected z map')
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
end


end