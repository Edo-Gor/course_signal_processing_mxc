%%
%   COURSE: Neural signal processing and analysis: Zero to hero
%  SESSION: Problem set: Spectral analyses of real and simulated data
%  TEACHER: Mike X Cohen, sincxpress.com
%

clear, clc

%% 1) Generate 10 seconds of data at 1 kHz, comprising 4 sine waves with different 
%    frequencies (between 1 and 30 Hz) and different amplitudes.
%    Plot the individual sine waves, each in its own plot. In separate subplots,
%    plot the summed sine waves with (1) a little bit of noise and (2) a lot of noise.
% 

srate  = 1000;
frex   = linspace(1,30,4);
amplit = [ 3 2 1 1 ];
phases = [ pi/4 pi/2 pi pi/3 ];
time   = 0:1/srate:10;

% create sine waves, first initialize to correct size
sine_waves = zeros(length(frex),length(time));

for fi=1:length(frex)
    sine_waves(fi,:) =  amplit(fi) * sin(2*pi*time*frex(fi) + phases(fi));
end

figure(1), clf

% plot constituent sine waves (without noise)
for snum=1:4
    subplot(4,1,snum)
    plot(time,sine_waves(snum,:))
    title([ 'Sine wave component with frequency ' num2str(frex(snum)) ' Hz and phase ' num2str(phases(snum)) ])
end
xlabel('Time (s)'), ylabel('Amplitude (arb.)')

% sum sinewaves and add noise
littleNoise = 0.6 * randn(1,length(time)) + sum(sine_waves,1);
lotsOfNoise = 10  * randn(1,length(time)) + sum(sine_waves,1);


% plot summed sine waves with little/lot of noise
figure(2), clf
subplot(211)
plot(time,littleNoise)
title('Time series with LITTLE noise')

subplot(212)
plot(time,lotsOfNoise)
title('Time series with A LOT of noise')

%% 2) Compute the power spectrum of the simulated time series (use FFT) and plot the results, 
%    separately for a little noise and a lot of noise. Show frequencies 0 to 35 Hz.
%    How well are the frequencies reconstructed, and does this depend on noise?
%    
%    Frequencies are clearly shown in both spectra, noise influences the
%    clearness of the graphs, but the frequency domain is clearly easier to
%    read than the time-domain

figure(3), clf

% normalisation factor
pnts = length(time);

for noisei=1:2
    
    % FFT
    if noisei==1
        f = 2*abs(fft(littleNoise)/pnts);
    else
        f = 2*abs(fft(lotsOfNoise)/pnts);
    end
    
    % compute frequencies in Hz
    hz = linspace(0,srate/2,floor(pnts/2)+1);
    
    % plot the amplitude spectrum
    subplot(2,1,noisei)
    plot(hz,f(1:length(hz)),'linew',1.5)
    xlabel('Frequencies (Hz)'), ylabel('Amplitude')
    set(gca,'xlim',[0 35],'ylim',[0 max(amplit)*1.2]) % is this ([0 1]) a good choice for x-axis limit?!?!?! nope
    
    if noisei==1
        title('FFT with LITTLE noise')
    else
        title('FFT with LOTS OF noise')
    end
end


%% 3) Compute the power spectrum of data from electrode 7 in the laminar V1 data. 
%    First, compute the power spectrum separately for each trial and then average the power 
%    results together. 
%    Next, average the trials together and then compute the power spectrum. 
%    Do the results look similar or different, and are you surprised? Why might they look 
%    similar or different?
% 
%    The averaging done after FFT on each trial produces more noises at
%    neighbouring frequencies across the entire spectum; albeit the
%    structure is generally preserved, the FFT on the trial average outputs
%    less power and a flatter spectrum, which nearly cancel out the ~50Hz
%    component visible in the first spectrum (i.e., most of the gamma
%    oscillations are not phase-locked in the signal, thus if you average
%    the signal, this component vanishes)

clear

% load the LFP data and increase precision
load v1_laminar.mat
csd = double(csd);

% pick which channel
chan2use = 7;
pnts = length(timevec); % = length(csd(1,:,1))

% FFT of all trials individually (note that you can do it in one line!) the
% dimension is 1 because of the squeeze
powspectSeparate = fft(squeeze(csd(chan2use,:,:)),[],1)/pnts;
% Then average the single-trial spectra together (average over trials, not over frequencies)
powspectSeparate = mean((2*abs(powspectSeparate)).^2,2); % trials --> dim. 2

% now average trials first, then take the FFT of the trial average
powspectAverage  = mean(csd(chan2use,:,:),3);
powspectAverage  = (2*abs(fft(powspectAverage)/pnts)).^2;

% frequencies in Hz
hz = linspace(0,srate/2,floor(length(timevec)/2)+1);


% now plot
figure(4), clf
set(gcf,'name',[ 'Results from electrode ' num2str(chan2use) ])
subplot(211)
plot(hz,powspectSeparate(1:length(hz)),'linew',1.5)
set(gca,'xlim',[0 100],'ylim',[0 30000])
xlabel('Frequency (Hz)'), ylabel('Power')
title(['Averaging done after FFT on each trial: electrode ' num2str(chan2use)]) 

subplot(212)
plot(hz,powspectAverage(1:length(hz)),'linew',1.5)
xlabel('Frequency (Hz)'), ylabel('Power')
set(gca,'xlim',[0 100],'ylim',[0 30000])
title(['FFT done on trial average: electrode ' num2str(chan2use)])

%% 4) Do the same as above but for electrode 1. 
%    How do these results compare to the results from channel 7, and does this depend on
%    whether you average first in the time-domain or average the individual power spectra?
%    ANATOMICAL NOTE: channel 7 is around L4; channel 1 is in the hippocampus.
% 
%    the components visible at ~50Hz for electrod 7 is clearly missing;
%    there is just a 1/f structure which is way sharper for the averaging-before-FFT
%    than for the single trial FFTs 

clear

% load the LFP data and increase precision
load v1_laminar.mat
csd = double(csd);

% pick which channel
chan2use = 1;
pnts = length(timevec); % = length(csd(1,:,1))

% FFT of all trials individually (note that you can do it in one line!) the
% dimension is 1 because of the squeeze
powspectSeparate = fft(squeeze(csd(chan2use,:,:)),[],1)/pnts;
% Then average the single-trial spectra together (average over trials, not over frequencies)
powspectSeparate = mean((2*abs(powspectSeparate)).^2,2); % trials --> dim. 2

% now average trials first, then take the FFT of the trial average
powspectAverage  = mean(csd(chan2use,:,:),3);
powspectAverage  = (2*abs(fft(powspectAverage)/pnts)).^2;

% frequencies in Hz
hz = linspace(0,srate/2,floor(length(timevec)/2)+1);


% now plot
figure(4), clf
set(gcf,'name',[ 'Results from electrode ' num2str(chan2use) ])
subplot(211)
plot(hz,powspectSeparate(1:length(hz)),'linew',1.5)
set(gca,'xlim',[0 100],'ylim',[0 95000])
xlabel('Frequency (Hz)'), ylabel('Power')
title(['Averaging done after FFT on each trial: electrode ' num2str(chan2use)]) 

subplot(212)
plot(hz,powspectAverage(1:length(hz)),'linew',1.5)
xlabel('Frequency (Hz)'), ylabel('Power')
set(gca,'xlim',[0 100],'ylim',[0 95000])
title(['FFT done on trial average: electrode ' num2str(chan2use)])

%% 5) Fourier transform from scratch!
% Hey, wouldn't it be fun to program the discrete-time Fourier transform
% from scratch! Yes, of course it would be. Let's do that.
% 
% Generate a 20-element vector of random numbers.
% Use the hints below to help you write the Fourier transform.
% Next, use the fft function on the same data to verify that your FT was accurate.

clear

N       = 20;          % length of sequence
signal  = randn(1,N);  % data
srate   = 100;
fTime   = (0:N-1)/N;   % "time" used in Fourier transform

% initialize Fourier output matrix
fourierCoefs = zeros(size(signal)); 

% loop over frequencies
for fi=1:N
    
    % create sine wave for this frequency
    fourierSine = exp( -1i*2*pi*(fi-1)*fTime );
    
    % compute dot product as sum of point-wise elements
    fourierCoefs(fi) = sum(fourierSine.*signal);
    
end

% divide by N to scale coefficients properly anc compute Hz frequencies
fourierCoefs = fourierCoefs / N;
hz = linspace(0,srate/2,floor(N/2)+1); 

figure(5), clf
subplot(211)
plot(signal)
xlim([1,20])
title('Data')

subplot(212)
plot(hz,abs(fourierCoefs(1:length(hz))*2),'*-')
xlabel('Frequency (Hz)')

% for comparison, use the fft function on the same data
fourierCoefsF = 2*abs(fft(signal)/N);

% plot the results on top. Do they look similar? (Should be identical!)
hold on
plot(hz,fourierCoefsF(1:length(hz)),'ro')
legend({'Manual FT';'FFT'})
title('Fourier trasform')

%% 6) zero-padding and interpolation
% Compute the power spectrum of channel 7 from the V1 dataset. Take the
% power spectrum of each trial and then average the power spectra together.
% But don't use a loop over trials! And use only the data from 0-.5 sec. 
% 
% What is the frequency resolution?
%   > native freq. resolution is ~2 Hz
% 
% Repeat this procedure, but zero-pad the data to increase frequency
% resolution. Try some different zero-padding numbers. At what multiple of the native nfft
% does the increased frequency resolution have no appreciable visible effect on the results?
%   > it makes the spectrum smoother
%   > the shape of the spectrum never changes a lot actually 

clear

% load the LFP data and increase precision
load v1_laminar.mat
csd = double(csd);

chan2use = 7;
pnts = length(timevec); % = length(csd(1,:,1))

tidx(1) = dsearchn(timevec',0);
tidx(2) = dsearchn(timevec',.5);

% set nfft to be multiples of the length of the data
nfft = 20 * (diff(tidx)+1);

powspect = fft(squeeze(csd(chan2use,(tidx(1):tidx(2)),:)),nfft,1)/pnts; % data --> dim. 1 because of squeeze
powspect = mean((2*abs(powspect)).^2,2); % trials --> dim. 2

hz = linspace(0,srate/2,floor(nfft/2)+1);

figure(6), clf
plot(hz,powspect(1:length(hz)),'k-o')
xlabel('Frequency (Hz)'), ylabel('Power')
set(gca,'xlim',[0 120])
title([ 'Frequency resolution is ' num2str(mean(diff(hz))) ' Hz' ])


%% 7) Poor man's filter via frequency-domain manipulations.
%  The goal of this exercise is to see how a basic frequency-domain filter
%  works: Take the FFT of a signal, zero-out some Fourier coefficients,
%  take take the IFFT. 
%  Here you will do this by generating a 1/f noise signal and adding 50 Hz
%  line noise.

clear

% Generate 1/f noise with 50 Hz line noise. 
srate = 1234;
npnts = srate*3;
time  = (0:npnts-1)/srate;

% the key parameter of pink noise is the exponential decay (ed)
ed = 50; % try different values!
as = rand(1,npnts) .* exp(-(0:npnts-1)/ed);
fc = as .* exp(1i*2*pi*rand(size(as)));

signal = real(ifft(fc)) * npnts;

% now add 50 Hz line noise
signal = signal + 1*sin(2*pi*50*time);

% compute its spectrum
hz = linspace(0,srate/2,floor(npnts/2)+1);
signalX = fft(signal);


%%% plot the signal and its amplitude spectrum
figure(7), clf
subplot(211)
plot(time,signal,'k')
xlabel('Time (s)'), ylabel('Activity')
title('Time domain')

subplot(212)
plot(hz,2*abs(signalX(1:length(hz)))/npnts,'k')
set(gca,'xlim',[0 80])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Frequency domain')

%% zero-out the 50 Hz component

% find the index into the frequencies vector at 50 hz
hz50idx = dsearchn(hz',50);

% create a copy of the frequency vector
signalXf = signalX; % f=filter

% zero out the 50 Hz component
signalXf(hz50idx) = 0;

% take the IFFT of the spectrum whose 50Hz components is removed
signalf = real(ifft(signalXf));

% take back FFT of filtered signal
signalXf = fft(signalf);


% plot on top of original signal
figure(8), clf
subplot(211), hold on
plot(time,signal,'k')
plot(time,signalf,'r')
xlabel('Time (s)'), ylabel('Activity')
title('Time domain')

subplot(212), hold on
plot(hz,2*abs(signalX(1:length(hz)))/npnts,'k')
plot(hz,2*abs(signalXf(1:length(hz)))/npnts,'ro-','markerfacecolor','r')
set(gca,'xlim',[0 80])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Frequency domain')


%%% QUESTION: Why do you need to take the real part of the ifft?
%   > because you forgot the negative frequencies you fool!
%   but there might also be some residual imag part due to computer
%   rounding errors even though theoretically they should all cancel out,
%   even in case you correctly removed both pos. and neg. frequecies)
% 
%%% QUESTION: Why doesn't this procedure get rid of the line noise?!?!?!?
%   > because you forgot the negative frequencies you fool!
% 


%% now fix the problem ;)

% Notice that the filter didn't work: It attenuated but did not eliminate
% the line noise. Why did this happen? Use plotting to confirm your
% hypothesis! Then fix the problem in this cell ;)

% find the index into the frequencies vector at 50 hz
hz50idx = dsearchn(hz',50);

% create a copy of the frequency vector
signalXf = signalX; % f=filter

% zero out the 50 Hz component
signalXf(hz50idx)       = 0;
signalXf(end-hz50idx+2) = 0; % hint: the negative freqs. (+2 to account for DC and indexing(zero freq = position one))

% take the IFFT
signalff = ifft(signalXf);

% take FFT of filtered signal
signalXf = fft(signalff);


% plot all three versions
figure(9), clf
subplot(211), hold on
plot(time,signal,'k')
plot(time,signalf,'b')
plot(time,signalff,'r')
xlabel('Time (s)'), ylabel('Activity')
legend({'Original';'Half-filtered';'Filtered'})
title('Time domain')

subplot(212), hold on
plot(hz,2*abs(signalX(1:length(hz)))/npnts,'k')
plot(hz,2*abs(signalXf(1:length(hz)))/npnts,'ro-','markerfacecolor','r')
set(gca,'xlim',[0 80])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Frequency domain')

%% 8) Exploring why the filter in #7 isn't a good filter
%   > basically because just zeroing out a frequency bin generates a sharp
%   nonstationarity in the signal, such as the one more evident in the
%   simulation here below

% number of time points
N = 1000;

% generate a flat Fourier spectra
fourspect1 = ones(N,1);

% copy it and zero-out some fraction
fourspect2 = fourspect1;
fourspect2(round(N*.1):round(N*.2)) = 0;

% create time-domain signals via IFFT of the spectra
signal1 = real(ifft(fourspect1));
signal2 = real(ifft(fourspect2));
time = linspace(0,1,N);


% and plot!
figure(10), clf
subplot(211)
plot(time,fourspect1, time,fourspect2,'linew',2)
set(gca,'ylim',[0 1.05])
xlabel('Frequency (a.u.)')
title('Frequency domain')
legend({'Flat spectrum';'With edges'})

subplot(212)
plot(time,signal1, time,signal2,'linew',2)
set(gca,'ylim',[-1 1]*.05)
xlabel('Time (a.u.)')
title('Time domain')
legend({'Flat spectrum';'With edges'})


%% done.
