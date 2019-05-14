function results = stats(data)
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

datat = data(:);
meanval = nanmean(datat);
if round(meanval,10) == 0 
    meanval = 0;
end
sdval = nanstd(datat);
medianval = nanmedian(datat);
minval = nanmin(datat);
maxval = nanmax(datat);
N = length(datat);
NaNCounts = sum(isnan(datat));

% results = [meanval medianval minval maxval N NaNCounts];
% results = table(minval, meanval, medianval, maxval, int64(N), int64(NaNCounts));
% results.Properties.VariableNames = {'Min', 'Mean', 'Median', 'Max', 'N', 'NaNs'};

results = struct();
results.min = minval;
results.meanval = meanval;
results.sd = sdval;
results.medianval = medianval;
results.maxval = maxval;
results.n = int64(N);
results.NaNCounts = int64(NaNCounts);

end