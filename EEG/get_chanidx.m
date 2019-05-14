function [chanidx] = get_chanidx(eeg,channames)
% get_chanidx 

if isstruct(eeg)
    chanidx = find(ismember(lower({eeg.chanlocs.labels}),lower(channames)));
elseif iscell(eeg)
    chanidx = find(ismember(lower(eeg),lower(channames)));
end
end

