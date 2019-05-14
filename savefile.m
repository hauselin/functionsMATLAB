function [] = savefile(outputDirectory,filename,varargin)
% outputDirectory: path/directory to save variables to
% filename: name of file
% varargin: variables to save (do not use quotation marks! see usage); options for saving
%
% Usage
% savefile('Desktop/Results','results.mat',variable1,variable2,variable3);
% Last modified by Hause Lin 19-12-18 23:20 hauselin@gmail.com

if ~exist(outputDirectory)
    mkdir(outputDirectory)
end

outputfullpath = fullfile(outputDirectory,filename);
disp(' ');
disp(['Saving to ' outputfullpath '...']);
opts = {};

for ind = 1:length(varargin)
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