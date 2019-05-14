function [channames] = get_channame(eeg,chanidx)
% get_channame 
channames = {eeg.chanlocs(chanidx).labels};
end

