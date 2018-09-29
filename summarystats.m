function stats = summarystats(data)
% summarystats computes descriptive statistics on input data
% stats = summarystats(data) returns descriptive statistics: mean,
% median, minimum, maximum, count (length), and number of NaNs
%
% INPUT: a 1D vector of numbers
%
% OUTPUT: a table containing min, mean, median, max, count, and
% number of NaNs
% 
% Code written by Hause Lin hauselin@gmail.com
% Last modified by Hause Lin 15-08-17 9:28 PM

%% Compute descriptive statistics
% data = varargin{1};

meanval = nanmean(data);
medianval = nanmedian(data);
minval = nanmin(data);
maxval = nanmax(data);
N = length(data);
NaNCounts = sum(isnan(data));

% stats = [meanval medianval minval maxval N NaNCounts];
stats = table(minval, meanval, medianval, maxval, int64(N), int64(NaNCounts));
stats.Properties.VariableNames = {'Min', 'Mean', 'Median', 'Max', 'N', 'NaNs'};

%% display output stats in order of output
% disp(['Mean is ' num2str(meanval) '.'])
% disp(['Median is ' num2str(medianval) '.'])
% disp(['Minimum is ' num2str(minval) '.'])
% disp(['Maximum is ' num2str(maxval) '.'])
% disp(['Number of numbers is ' num2str(N) '.'])
% disp(['Number of NaNs is ' num2str(NaNCounts) '.'])

end