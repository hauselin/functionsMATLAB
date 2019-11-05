function [] = cnt_to_mat(cfg)
% cnt_to_mat
% converts EEProbe cnt data to eeglab EEG structure and saves as mat file

if nargin == 0 % if no input provided, return function argument information
    disp('cnt_to_mat function parameters: ');
    struct('rawdatadir','raw data path/directory',...
        'matdatadir','output directory path for mat files')
    return
end

if ~exist(cfg.matdatadir)
    mkdir(cfg.matdatadir);
end

files = dir(fullfile(cfg.rawdatadir,'*.cnt'));
    
for f=1:length(files)
    outputfile = fullfile(cfg.matdatadir,strrep(files(f).name,'.cnt','.mat'));
    if ~exist(outputfile)
        disp(['Converting ' files(f).name]);
        EEG = pop_loadeep_v4(fullfile(cfg.rawdatadir,files(f).name));
        save(outputfile,'EEG');
        clear EEG
    else
        disp(['Output file already exists: ' outputfile]);
    end
end

end

