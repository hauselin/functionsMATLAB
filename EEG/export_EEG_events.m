function [] = export_EEG_events(EEG)

% save triggers/event info
EEG.EVENTLIST = [];
EEG = pop_creabasiceventlist(EEG,'AlphanumericCleaning','on','BoundaryNumeric',{-99},'BoundaryString',{'boundary'});
el = EEG.EVENTLIST.eventinfo; 
elT = struct2table(el);
elT = elT(:,{'item','code','spoint'});
C = cell(size(elT,1),1); % empty cell to store subject id
C(:) = {EEG.subject}; % fill empty cell with subject id
elT = [C elT]; % join cell with table
elT.Properties.VariableNames = {'subject','eventN','eventCode','samplingPoint'}; % rename variable names

outDir = 'DesignMatrix';
if ~exist(outDir)
    mkdir(outDir);
end
writetable(elT,fullfile(outDir,[EEG.subject '_events.csv'])) % all events

end

