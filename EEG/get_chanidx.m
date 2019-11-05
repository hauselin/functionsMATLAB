function [chanidx] = get_chanidx(eeg,channames)
% get_chanidx(eeg,channames) 
% returns channel indices/numbers in eeg (EEGLAB structure) for channels in channames
% 
% eeg: EEG data structure (requires eeg.chanlocs field)
% channames: channel names to get indices of

if isstruct(eeg)
    chanidx = find(ismember(lower({eeg.chanlocs.labels}),lower(channames)));
elseif iscell(eeg)
    chanidx = find(ismember(lower(eeg),lower(channames)));
end
end

