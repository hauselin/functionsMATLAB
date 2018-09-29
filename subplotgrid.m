function rowcol = subplotgrid(totalSubplots, currentSubplot)
% subplotgrid determines and returns rows and columns needed for a given number of plots (totalSubplots).
% Plots the subplot at position currentSubplot
%
% USAGE
%   rowcol = subplotgrid(subplots, currentSubplot) % returns 1 x 2 vector
%   rowcol = subplotgrid(5, 3) % 5 subplots in total; make the 3rd plot active
% 
% INPUTs: 
%   totalSubplots: integer indicating total number of subplots
%   currentSubplot: integer indicating the current active plot    
%
% OUTPUTS
%   rowcol: 1 x 2 vector containing row and column numbers to draw
% 
% Code written by Hause Lin hauselin@gmail.com 
% Last modified by Hause Lin 20-09-17 17:03

%% Determine row and columns to draw on figure

r = ceil(totalSubplots /ceil(sqrt(totalSubplots))); % no. of rows
c = ceil(sqrt(totalSubplots)); % no. of cols
rowcol = [r c];

subplot(r, c, currentSubplot) % plots figure on the correct subplot
    
end