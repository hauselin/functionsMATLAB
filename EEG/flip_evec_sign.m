function [activationpatterns, timeseriescomponents, evecsignflip] = flip_evec_sign(evecs,chan,activationpatterns,timeseriescomponents,chanlocs)
% flip_evec_sign
% Flips sign of eigenvectors from GED/PCA based on the sign at one EEG channel
%
% INPUTS
% evecs: EEGchan by eigenvector matrix
% chan: which chan to use for sign flipping (e.g., 'FCz')
% activationpatterns: EEGchan by activation patterns from GED 
% timeseriescomponents: eigenvector by time matrix (component time series)
% chanlocs: eeglab structure chanlocs field
%
% Last modified by Hause Lin 15-05-19 0:10 PM hauselin@gmail.com

chanidx = find(strcmpi({chanlocs.labels},chan));
evecsignflip = sign(activationpatterns(chanidx,:)); % sign at channel for each activation pattern/topography

% flip activation pattern based on sign at channel
activationpatterns = repmat(evecsignflip,size(activationpatterns,1),1) .* activationpatterns;

% flip time series based on sign at channel
ncomps = size(timeseriescomponents,1);
timepoints = size(timeseriescomponents,2);

timeseriescomponents = timeseriescomponents .* repmat(evecsignflip(1:ncomps)',1,timepoints);
end

