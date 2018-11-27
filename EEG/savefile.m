function [] = savefile(outputDirectory,filename,varargin)

if ~exist(outputDirectory)
    mkdir(outputDirectory)
end

outputfullpath = fullfile(outputDirectory,filename);
disp(' ');
disp(['Saving to ' outputfullpath '...']);
opts = {};

for ind = 1:length(varargin)
    % ind = 1
    currentVariableName = inputname(ind+2);
    % Handle save options such as -append, -v7.3, etc
    if strcmp(currentVariableName,'')
        opts{length(opts) + 1} = varargin{ind};
    else
        varNames.(currentVariableName) = varargin{ind};
        disp(['Saving ' currentVariableName '...']);
    end
end

if ~isempty(opts)
    save(outputfullpath,'-struct','varNames',opts{:})
else
    save(outputfullpath,'-struct','varNames')
end

end