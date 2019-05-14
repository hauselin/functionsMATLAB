function [EEG] = rerun_iclabel(EEG)

eeglab redraw;
EEG = iclabel(EEG);
nbchan = min(size(EEG.icawinv));
if nbchan > 35
    nbchan = 30;
end
pop_viewprops(EEG,0,1:nbchan,{'freqrange',[2 70]},{},1,'ICLabel');

end

