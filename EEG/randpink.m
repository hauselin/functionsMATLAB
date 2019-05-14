function [noiseout] = randpink(n,varargin)
% Usage
% randpink(100,50,'plot') % generate 10 values, exponential decay 50, plot
% randpink(100) % generate 10 values, exponential decay 20, no plot
% Last modified by Hause Lin 03-04-19 1:28 PM hauselin@gmail.com

switch nargin
    case 1 % no options
        ed = 20;
        showplot = false;
    case 2
        ed = varargin{1};
        showplot = false;
    case 3
        ed = varargin{1};
        showplot = varargin{2};
        if strcmpi(showplot,'plot')
            showplot = true;
        else
            showplot = false;
        end
end

srate = 500; % sampling rate in Hz
time = linspace(0,srate,n+1)/(srate);
pnts = length(time);
hz = linspace(0,srate/2,floor(length(time)/2)+1);

% generate 1/f amplitude spectrum
as = rand(1,floor(pnts/2)-1) .* exp(-(1:floor(pnts/2)-1)/ed); % scaling noise by 1/(f^ed) (pointwise multiply)
as = [as(1) as 0 0 as(:,end:-1:1)]; % mirror image

% Fourier coefficients
fc = as .* exp(1i*2*pi*rand(size(as))); % add random phases

% inverse Fourier transform to create the noise
noise = real(ifft(fc)) * pnts;
noiseout = noise(1:n);

if showplot
    figure(100), clf
    subplot(211)
    plot(time,noise,'k')
    set(gca,'fontsize',15)
    title('Pink noise: Time domain')
    xlabel('Time (s)'), ylabel('Amplitude')

    subplot(223)
    [y,x] = hist(noise,100);
    plot(x,y,'k','linew',2)
    xlabel('Values'), ylabel('N per bin')
    title('Signal histogram (distribution)')
    set(gca,'fontsize',15)

    subplot(224)
    amp = abs(fft(noise)/pnts);
    amp(2:end) = 2*amp(2:end);
    plot(hz,amp(1:length(hz)),'k')
    title('Frequency domain')
    set(gca,'fontsize',15)
    xlabel('Frequency (Hz)'), ylabel('Amplitude')
end

end