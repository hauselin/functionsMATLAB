function [indicesStruct] = get_epochidx(cfg)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%% check consistency

if cfg.data.trials ~= size(cfg.designmat,1)
    error('Number of epochs in data and design matrix do not match!');
end

%% get column indices

indices = [];
indicesStruct = struct();
indicesStruct.idx = indices;
indicesStruct.conds = cfg.filter;

% filter table to get row indices
for f=1:length(cfg.filter)
    if ischar(cfg.filter{f}) % if filter is character
        matchingindices = find(strcmp(cfg.designmat.(cfg.var),cfg.filter{f})); % find matching strings
    elseif isnumeric(cfg.filter{f}) % filter is numeric
        matchingindices = find(cfg.designmat{:,cfg.var} == cfg.filter{f}); % find matching numbers
    end
    indicesStruct.(['c',num2str(f)]) = matchingindices;
    indices = [indices; matchingindices];
end

indicesStruct.idx = indices;
indicesStruct.all = [cfg.designmat(indices,[cfg.var]) num2cell(indices)];
indicesStruct.designmat = cfg.designmat(indices,:);

%% print results

disp(['Epochs matched: ' num2str(length(indices))]);

end

