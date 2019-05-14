% reorderchans
% Reorders EEG channel/sensor.
%
% Usage examples
% new = {'Fp1','Fpz','Fp2','F8','F4','Fz','F3','F7','FC3','FCz','FC4','C4','Cz','C3','FT7','T3','FT8','T4','CP4','CPz','CP3','P3','Pz','P4','TP8','T6','TP7','T5','O1','Oz','O2','Corr','Zygo'};
% EEGsorted = reorderchans(EEG,new,true);
%
% Arguments
% reorderchans(EEG,neworder,showplot) 
% EEG: EEG structure from eeglab
% neworder: new order of channels (cell array)
% showplot: whether to show figures (true or false)
%
% Written in MATLAB 2017b. Last modified by Hause Lin 03-08-18 08:30 hauselin@gmail.com

function [EEGout] = reorderchans(EEG,neworder,showplot)

EEGout = EEG; % make output EEG structure
for chanI = 1:length(EEGout.chanlocs) % assign new indices (1 to n)
    EEGout.chanlocs(chanI).urchan = chanI;
end

% check if neworder length and number of channels in EEG match
if (EEG.nbchan ~= length(neworder))
    error('Number of channels in EEG and neworder should be the same!');
end

% add channel information if necessary
if sum([EEGout.chanlocs.X]) == 0 % no channel info
    EEGout = pop_chanedit(EEGout,'load',[],'lookup','standard-10-5-cap385.elp'); % add channel info
end

if showplot && sum([EEGout.chanlocs.X]) ~= 0 
    figure(1029482)
    clf
    set(gcf,'name','Channel order','numbertitle','off')
    subplot(241)
    topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'numpoint', 'chaninfo', EEG.chaninfo);
    title('Original location numbers')
    subplot(242)
    topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    title('Original location labels')
    
    subplot(243)
    if ndims(EEG.data) == 2
        imagesc(squeeze(nanmean(EEG.data(:,1:10),2)))
    elseif ndims(EEG.data) == 3
        imagesc(squeeze(nanmean(EEG.data(:,1:10,1),2)))
    end
    
    subplot(244)
    if ndims(EEG.data) == 2
        topoplotIndie(squeeze(nanmean(EEG.data(:,1:10),2)),EEG.chanlocs,'electrodes','on');
    elseif ndims(EEG.data) == 3
        topoplotIndie(squeeze(nanmean(EEG.data(:,1:10,1),2)),EEG.chanlocs,'electrodes','on');
    end
        
    title('Original chan mean activity')
    colorbar
    colormap parula
%     for a = 1:length({EEG.chanlocs.labels})
%         outlabel{a} = strcat(EEG.chanlocs(a).labels,['-' num2str(a)]);
%     end
    set(gca,'ytick',1:EEG.nbchan,'yticklabel',{EEG.chanlocs.labels},'xtick',1,'xticklabel','Mean amplitude (first 10 points)')
end

for chanI = 1:length(EEGout.chanlocs) % add neworder field in EEG.chanlocs
    newPosition = find(strcmpi(neworder,EEGout.chanlocs(chanI).labels));
    EEGout.chanlocs(chanI).neworder = newPosition;
end

% fill in empty field values with NaN
for i = 1:size(EEGout.chanlocs,2)
    if isempty(EEGout.chanlocs(i).neworder)
        EEGout.chanlocs(i).neworder = NaN;
    end
end

[sorted, order] = sort([EEGout.chanlocs.neworder],'ascend'); % sort cell array
sortedChans = EEGout.chanlocs(order);
sortChanIdx = [sortedChans.urchan]; % indices to use to resort EEG.data later on (MOST IMPORTANT THING!)

for chanI = 1:length(sortedChans) % assign new indices (1 to n)
    sortedChans(chanI).neworder = chanI;
    sortedChans(chanI).urchan = chanI; % replace urchan with new chan index
end

EEGout.chanlocs = sortedChans; % replace with sorted chanlocs

% resort data
if ndims(EEG.data) == 2
    EEGout.data = EEG.data(sortChanIdx,:); % resort data
elseif ndims(EEG.data) == 3
    EEGout.data = EEG.data(sortChanIdx,:,:); % resort data
else
    error('check EEG.data dimensions!');
end 

if showplot && sum([EEGout.chanlocs.X]) ~= 0 
    subplot(245)
    topoplot([],EEGout.chanlocs, 'style', 'blank',  'electrodes', 'numpoint', 'chaninfo', EEG.chaninfo);
    title('Sorted location numbers')
    subplot(246)
    topoplot([],EEGout.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    title('Sorted location labels')
    
    subplot(247)
    if ndims(EEG.data) == 2
        imagesc(squeeze(nanmean(EEGout.data(:,1:10),2)))
    elseif ndims(EEG.data) == 3
        imagesc(squeeze(nanmean(EEGout.data(:,1:10,1),2)))
    end
    
    subplot(248)
    if ndims(EEG.data) == 2
        topoplotIndie(squeeze(nanmean(EEGout.data(:,1:10),2)),EEGout.chanlocs,'electrodes','on');
    elseif ndims(EEG.data) == 3
        topoplotIndie(squeeze(nanmean(EEGout.data(:,1:10,1),2)),EEGout.chanlocs,'electrodes','on');
    end
    
    title('Sorted chan mean activity')
    colorbar
    try
        colormap viridis
    catch
        colormap parula
    end
            
    set(gca,'ytick',1:EEG.nbchan,'yticklabel',{EEGout.chanlocs.labels},'xtick',1,'xticklabel','Mean amplitude (first 10 points)')
    set(findall(gcf,'-property','FontSize'),'FontSize',11)
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'PaperPositionMode','auto','DefaultTextInterpreter','none','PaperOrientation','portrait'); % maximize figure
end

% return cell for checking purposes
chanActivityOriginal = cell(2,length({EEGout.chanlocs.labels}));
chanActivitySorted = chanActivityOriginal;

chanActivityOriginal(1,:) = {EEG.chanlocs.labels};
chanActivitySorted(1,:) = {EEGout.chanlocs.labels};

if ndims(EEG.data) == 2
    d1 = round(nanmean(EEG.data(:,1:10,:),2)',3);
    d2 = round(nanmean(EEGout.data(:,1:10,:),2)',3);
elseif ndims(EEG.data) == 3
    d1 = round(nanmean(EEG.data(:,1:10,1),2)',3);
    d2 = round(nanmean(EEGout.data(:,1:10,1),2)',3);        
end

for d = 1:length(d1)
    chanActivityOriginal{2,d} = sprintf('%.2f',d1(d));
    chanActivitySorted{2,d} = sprintf('%.2f',d2(d));
end

disp('Mean activity for first 10 points in each channel');
disp('Cell array columns');
disp('1. original order; 2. original activity; 3. requested order; 4. sorted order; 5. sorted activity');
[chanActivityOriginal' neworder' chanActivitySorted']

% sortrows([chanActivityOriginal' chanActivitySorted'],3)

end