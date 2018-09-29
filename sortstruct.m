function[structout] = sortstruct(structure,fieldname,order)

% fill in empty field values with NaN
for i = 1:size(structure,2)
    if isempty(structure(i).(fieldname))
        structure(i).(fieldname) = NaN;
    end
end

if isa(structure(1).(fieldname),'char') 
    cellToSort = {structure.(fieldname)};
    [sorted, order] = sortrows(cellToSort',order); % sort cell array
    structout = structure(order);
else
    cellToSort = {structure.(fieldname)};
    [sorted, order] = sort(cellToSort,order); % sort cell array
    structout = structure(order);

% remove NaNs
for i = 1:size(structout,2)
    if isnan(structout(i).(fieldname))
        structout(i).(fieldname) = [];
    end
end

end
