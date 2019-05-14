function [results] = gedSliding_eeg(cfg)
% gedSliding_eeg
% Performs GED on a sliding window. Calls ged_eeg function to perform each GED.
%
% Type the function name without any inputs to get parameter information
% 
% Last modified by Hause Lin 27-04-19 8:59 PM hauselin@gmail.com

if nargin == 0 % if no input provided, return function argument information
    disp('ged_eeg function parameters: ');
    struct('data','eeglab data structure',...
        'singletrial','single-trial covariance (1 or 0; default: 1, yes)',...
        'Sstarttime','S matrix sliding window start time',...
        'Swinsize','S matrix window size',...
        'Sstepsize','S matrix sliding window step size',...
        'Rwin','R matrix time window',...
        'componentskeep','number of components/maps to keep',...
        'verbose','print messages for debugging (default: 0)')
    return
end

if ~isfield(cfg,'singletrial')
    cfg.singletrial = 1;
end

if ~isfield(cfg,'Sendtime')
    cfg.Sendtime = max(cfg.data.times);
end

if ~isfield(cfg,'verbose')
    cfg.verbose = 0;
end

%% define sliding time windows

% S matrices
Sslidewindows = cfg.Sstarttime:cfg.Sstepsize:(cfg.Sendtime-cfg.Swinsize);
Sslidewindows = [Sslidewindows' (Sslidewindows + cfg.Swinsize)'];

% R matrix is fixed

%% define output matrices

slide_eigvals = zeros(cfg.data.nbchan,length(Sslidewindows)); % eigenval_win
slide_eigvalsprop = slide_eigvals;
slide_eigvecs = zeros(cfg.data.nbchan,cfg.data.nbchan,length(Sslidewindows)); % eigenvec_eigenvec_win
slide_covS = slide_eigvecs; % chan_chan_win
slide_covR = zeros(cfg.data.nbchan); % chan_chan_win
slide_activationPattern = zeros(cfg.componentskeep,cfg.data.nbchan,length(Sslidewindows)); % activation_component_win

%% loop and call ged_eeg() function

disp(['GED on ' num2str(size(Sslidewindows,1)) ' windows...']);

for wini=1:size(Sslidewindows,1) % for each window
    if cfg.verbose
        disp(['S window ' num2str(wini) ': ' num2str(Sslidewindows(wini,1)) ' to ' num2str(Sslidewindows(wini,2))]);
    end
    tempcfg = []; 
    tempcfg.singletrial = cfg.singletrial; 
    tempcfg.data = cfg.data;
    tempcfg.Rwin = [cfg.Rwin(1) cfg.Rwin(2)]; 
    tempcfg.Swin = [Sslidewindows(wini,1) Sslidewindows(wini,2)];
    tempcfg.activationcomponents = cfg.componentskeep; 
    tempcfg.timeseriescomponents = cfg.componentskeep;
    tempcfg.verbose = cfg.verbose;
    tempres = ged_eeg(tempcfg);
    
    slide_eigvals(:,wini) = tempres.evals;
    slide_eigvalsprop(:,wini) = tempres.evalsprop;
    slide_eigvecs(:,:,wini) = tempres.evecs;
    slide_covS(:,:,wini) = tempres.S;
    slide_covR = tempres.R;
    slide_activationPattern(:,:,wini) = tempres.activationpatterns';
end

disp(['Done!']);

%% return result

results = [];
results.Swins = Sslidewindows;
results.S = slide_covS;
results.R = slide_covR;
results.evecs = slide_eigvecs;
results.evals = slide_eigvals;
results.evalsprop = slide_eigvalsprop;
results.activationpatterns = slide_activationPattern;
results.chanlocs = cfg.data.chanlocs;

end

