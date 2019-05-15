function [peakinfo] = get_peakinfo(data,times,timewindow,minmax,border,plotfigure)
% get_peakinfo 
% Locates peak (min or max) within a time window (only works for 1D input)
% 
% USAGE
% get_peakinfo(EEG.data,EEG.times,[150 350],'max',20,true)
% 
% INPUTS
% data: structure that contains the data field (e.g., EEG.data)
% times: time points corresponding to time dimension of input data
% timewindow: start and end times to find peak
% minmax: 'min' or 'max' value to identify
% border: how much time to "pad" around the peak 
% plotfigure: whether to plot results
%
% Last modified by Hause Lin 11-05-19 11:59 AM hauselin@gmail.com

% check if lengths match
if ~length(data) == length(times)
    error('Length of data and times have to match!');
end

peakinfo = struct();

timewinidx = dsearchn(times',timewindow');
dataNaN = data;
dataNaN(1:timewinidx(1)-1) = NaN;
dataNaN(timewinidx(2):end) = NaN;

% get max info
[maxval maxidx] = max(dataNaN);

% get min info
[minval minidx] = min(dataNaN);

switch minmax
    case 'max'
        peakinfo.peakval = maxval;
        peakinfo.peaktime = times(maxidx);
        peakinfo.peakindex = maxidx;
    case 'min'
        peakinfo.peakval = minval;
        peakinfo.peaktime = times(minidx);
        peakinfo.peakindex = minidx;
end

if plotfigure
    figure
    subplot(211)
    plot(times,data,'linew',1);
    hold on
    plot(times,dataNaN,'linew',2);
    y = ylim;
    plot([peakinfo.peaktime peakinfo.peaktime],[y(1) y(2)],'linew',2);
    if border > 0
        rectangle('Position',[peakinfo.peaktime-border y(1) border*2 diff(y)],'FaceColor',[0.5 0.5 0.5 0.4],'linestyle','none');
    end
    title([minmax ': ' num2str(peakinfo.peakval) ', ' minmax ' time: ' num2str(peakinfo.peaktime)]);
    
    subplot(212)
    plot(times,data,'linew',1);
    hold on
    plot(times,dataNaN,'linew',2);
    y = ylim;
    plot([peakinfo.peaktime peakinfo.peaktime],[y(1) y(2)],'linew',2);
    if border > 0
        rectangle('Position',[peakinfo.peaktime-border y(1) border*2 diff(y)],'FaceColor',[0.5 0.5 0.5 0.4],'linestyle','none');
    end
    xlim([timewindow(1)-border*2 timewindow(2)+border*2]);
end

peakinfo.peakstart = peakinfo.peaktime-border;
peakinfo.peakend = peakinfo.peaktime+border;
peakinfo.duration = peakinfo.peakend - peakinfo.peakstart;

end

