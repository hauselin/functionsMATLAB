function [result] = matchvectors(template,vectors,measure)
%matchvectors Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    measure = "all";
end

% simple correlations
rs = corrcoef([template vectors]);
rs = rs(2:end,1)';

% dot products
dps = template' * vectors;

% euclidean distance
eds = sum((repmat(template,1,size(vectors,2)) - vectors).^2,1);

% cosine distance
cosds = pdist([template'; vectors'],'cosine');
cosds = 1-cosds(1:size(vectors,2));

eds = sum((repmat(template,1,size(vectors,2)) - vectors).^2,1);

% all similarity indices
all = [rs; dps; -eds; cosds]; % reverse euclidean distance so larger = better

% get max value
[maxval,maxidx] = max(all,[],2);

% get modal values
% maxidx = [2 2 1 1 1 3 3 3]
% maxidx = [1 2 1 1 3 3 1 3 1]
U = unique(maxidx);
H = histc(maxidx, U);
Modes = U(H == max(H));

disp('Rows: correlation, dot product, euclidean, cosine');
disp([maxidx all]);

if strcmpi(measure,'all')
    if length(Modes) > 1 % if more than one unique value (different measures disagree)
        result = maxidx(2); %
        disp(['More than one mode! Returning best index based on dotproduct: ' num2str(result)]);
    else % all measures return same value (just 1 mode)
        result = mode(maxidx);
        disp(['Returning index of modal max value: ' num2str(result)]);
    end
elseif strcmpi(measure,'correlation')
    result = maxidx(1);
elseif strcmpi(measure,'dotproduct')
    result = maxidx(2); 
elseif strcmpi(measure,'euclidean')
    result = maxidx(3); 
elseif strcmpi(measure,'cosine')
    result = maxidx(4);
end

end

