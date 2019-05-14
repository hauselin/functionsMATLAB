function [] = save_figure(outputDirectory,filename,varargin)

if ~exist(outputDirectory)
    mkdir(outputDirectory)
end

outputfullpath = fullfile(outputDirectory,filename);
disp(' ');
disp(['Saving to ' outputfullpath '...']);

% save figure
print(gcf,'-djpeg','-r300', outputfullpath);

end