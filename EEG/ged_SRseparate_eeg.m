function [results] = ged_SRseparate_eeg(cfg)
% ged_SRseparate_eeg 
% Performs GED with S and R matrices from separate datasets
%
% Type the function name without any inputs to get parameter information
%
% Last modified by Hause Lin 27-04-19 8:59 PM hauselin@gmail.com

if nargin == 0 % if no input provided, return function argument information
    disp('ged_SRseparate_eeg function parameters: ');
    struct('Sdata','eeglab data structure (e.g., filtered data)',...
        'Rdata','eeglab data structure (e.g., broadband data)',...
        'Swin','S matrix time window',...
        'Rwin','R matrix time window',...
        'singletrial','single-trial covariance (1 or 0; default: 1, yes)',...
        'activationcomponents','number of activation components to plot',...
        'timeseriescomponents','number of time series components to plot',...
        'timeseriescomponents_with_SRdata','use Sdata or Rdata to compute time series components (default: "Rdata")',...
        'plotfig','plot results figure number (default: 0)',...
        'plottime','plot xlim for time series (default: [min max])',...
        'verbose','print messages for debugging (default: 0)')
    return
end

if ~isfield(cfg,'singletrial')
    cfg.singletrial = 1;
end

if ~isfield(cfg,'plotfig')
    cfg.plotfig = 0;
end

if ~isfield(cfg,'verbose')
    cfg.verbose = 0;
end

if ~isfield(cfg,'timeseriescomponents_with_SRdata')
    cfg.timeseriescomponents_with_SRdata = 'Rdata';
end

%% checks

if sum(size(cfg.Sdata.data) - size(cfg.Rdata.data)) ~= 0
    error('Sdata and Rdata dimensions do not match!');
end

%% get time indices for R and S matrices

cfg.Rwinidx = dsearchn(cfg.Rdata.times',[cfg.Rwin(1) cfg.Rwin(end)]');
cfg.Swinidx = dsearchn(cfg.Sdata.times',[cfg.Swin(1) cfg.Swin(end)]');

%% perform GED

cfg.Sdata.data = double(cfg.Sdata.data); % ensure double precision for eigendecomposition
cfg.Rdata.data = double(cfg.Rdata.data); % ensure double precision for eigendecomposition

if cfg.verbose
    disp(['Performing GED on ' num2str(cfg.Sdata.trials) ' trials']);
end

[covR,covS] = deal(zeros(cfg.Sdata.nbchan)); % initialize empty covariance matrices

if cfg.singletrial % 
    if cfg.verbose
        disp('Computing single-trial covariance matrices...');
    end
    for ti=1:cfg.Sdata.trials % covariance matrix for each trial
        % R matrix
        tdat = squeeze(cfg.Rdata.data(:,cfg.Rwinidx(1):cfg.Rwinidx(2),ti)); 
        covR = covR + cov(tdat');

        % S matrix
        tdat = squeeze(cfg.Sdata.data(:,cfg.Swinidx(1):cfg.Swinidx(2),ti)); 
        covS = covS + cov(tdat');
    end
    covR = covR./cfg.Rdata.trials; % average covariance
    covS = covS./cfg.Sdata.trials;
else % covariance matrix on avg data
    if cfg.verbose
        disp('Computing covariance matrices on averaged data...');
    end
    avgdataR = mean(cfg.Rdata.data,3);
    avgdataS = mean(cfg.Sdata.data,3);
    % R matrix
    tdat = squeeze(avgdataR(:,cfg.Rwinidx(1):cfg.Rwinidx(2))); 
    covR = cov(tdat');

    % S matrix
    tdat = squeeze(avgdataS(:,cfg.Swinidx(1):cfg.Swinidx(2))); 
    covS = cov(tdat');
end

% sort eigenvectors
[evecs,evals] = eig(covS,covR);
if ~isreal(evals)
    error('eigendecomposition returned imaginary values!');
end
[evals,sidx] = sort(diag(evals),'descend');
evecs = evecs(:,sidx);
evalsprop = evals ./ sum(evals) * 100;

%% compute filter-forward model/activation pattern and flip sign

% compute activation pattern/topography (a = wS, where S is covariance matrix S)
activationpatterns = [];
if cfg.activationcomponents
    activationpatterns = zeros(cfg.Sdata.nbchan,cfg.activationcomponents);
    for c=1:cfg.activationcomponents
        % compute activation pattern/topography (a = wS, where S is covariance matrix S)
        activationpatterns(:,c) = evecs(:,c)'*covS; % get component
        [~,idx] = max(abs(activationpatterns(:,c)));  % find max magnitude
        activationpatterns(:,c) = activationpatterns(:,c)*sign(activationpatterns(idx,c)); % possible sign flip
    end    
end

%% compute filter forward model and flip sign

evecs = evecs./sqrt(sum(evecs.^2,1));
% compute component time series (projections) (c = wX, where X is data matrix)
timeseriescomponents = [];
if cfg.timeseriescomponents
    if strcmpi(cfg.timeseriescomponents_with_SRdata,'Rdata')
        timeseriescomponents = evecs(:,1:cfg.timeseriescomponents)'*squeeze(mean(cfg.Rdata.data,3)); % comp_chan * chan_time
    elseif strcmpi(cfg.timeseriescomponents_with_SRdata,'Sdata')
        timeseriescomponents = evecs(:,1:cfg.timeseriescomponents)'*squeeze(mean(cfg.Sdata.data,3)); % comp_chan * chan_time
    end    
end

%% visualize results

if cfg.plotfig
    
    figure(cfg.plotfig)
    
    % plot eigenspectrum
    subplot(3,cfg.activationcomponents,1:(cfg.activationcomponents-2))
    plot(1:length(evalsprop),evalsprop,'s-','linew',1,'markersize',5,'markerfacecolor','k')
    set(gca,'xlim',[1 length(evalsprop)])
    ylabel('Variance explained'); title(['Eigenvalues']);
    if length(evalsprop) > 30
        title(['Eigenvalues for top 30 (of ' num2str(length(evalsprop)) ' values)']);
        set(gca,'xlim',[1 30])
    end
    
    subplot(3,cfg.activationcomponents,cfg.activationcomponents-1)
    tempdat = cfg.Sdata.data(:,cfg.Swinidx(1):cfg.Swinidx(2),:);
    topoplotIndie(squeeze(mean(mean(tempdat,2),3)),cfg.Sdata.chanlocs);
    title('Data for GED S matrix');
    
    subplot(3,cfg.activationcomponents,cfg.activationcomponents)
    tempdat = cfg.Rdata.data(:,cfg.Rwinidx(1):cfg.Rwinidx(2),:);
    topoplotIndie(squeeze(mean(mean(tempdat,2),3)),cfg.Rdata.chanlocs);
    title('Data for GED R matrix');
   
    % plot activation patterns/topography (a = wS)
    for i=1:cfg.activationcomponents
        subplot(3,cfg.activationcomponents,cfg.activationcomponents+i)
        topoplotIndie(activationpatterns(:,i),cfg.Sdata.chanlocs,'numcontour',0,'electrodes','off');
        title([ 'Component ' num2str(i) ])
        % axis tight
        colorbar
        caxis([-max(abs(caxis)) max(abs(caxis))])
    end
    
    % plot component time series (projections) (c = wX, where X is data matrix)
    for i=1:cfg.activationcomponents
        subplot(3,cfg.activationcomponents,2*cfg.activationcomponents+i)
        plot(cfg.Sdata.times,timeseriescomponents(i,:));
        if isfield(cfg,'plottime')
            xlim([cfg.plottime(1) cfg.plottime(2)])
        else
            xlim([cfg.Sdata.times(1) cfg.Sdata.times(end)])
        end
        xlabel('Time (s)')
        if i == 1 
            ylabel('Amplitude');
        end
        % axis tight
    end
    try
        colormap viridis
    end
end

%% return result

results = [];
results.R = covR;
results.S = covS;
results.evecs = evecs;
results.evals = evals;
results.evalsprop = evalsprop;
results.activationpatterns = activationpatterns;
results.timeseriescomponents = timeseriescomponents;
results.times = cfg.Sdata.times;
results.chanlocs = cfg.Sdata.chanlocs;

end

