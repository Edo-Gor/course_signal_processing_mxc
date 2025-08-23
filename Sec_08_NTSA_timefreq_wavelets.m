%%
%     COURSE: Solved problems in neural time series analysis
%    SECTION: Time-frequency analyses
%      VIDEO: Create a family of complex Morlet wavelets
% Instructor: sincxpress.com
%

clear, clc

%%

% time parameters
srate = 1025;
npnts = 2001; % use an odd number to center around 0!
timev = linspace(-1,1,npnts) * (npnts-1)/srate/2; % scale by sampling rate

% wavelet parameters
minFreq =  2; % Hz
maxFreq = 54; % Hz
numfrex = 99;

% vector of wavelet frequencies, log spaced
frex = logspace(log10(minFreq),log10(maxFreq),numfrex);

% Gaussian parameters
numCycl = logspace(log10(3),log10(15),numfrex);
fwhms   = logspace(log10(1),log10(.3),numfrex);

%% now create the wavelets

wavefam = zeros(2,numfrex,npnts);

for fi=1:numfrex
    
    % create complex sine wave
    csw = exp( 1i* 2*pi*frex(fi)*timev );
    
    % create the two Gaussians, one with cycles formula, one with FWHM formula
    s = numCycl(fi) / (2*pi*frex(fi));  % normalise cycles
    gaus1 = exp( -timev.^2 / (2*s)^2 );
    
    gaus2 = exp( -4*log(2)*timev.^2 / fwhms(fi)^2 );
    
    % now create the complex Morlet wavelets
    wavefam(1,fi,:) = csw .* gaus1;
    wavefam(2,fi,:) = csw .* gaus2;
end

%% image the wavelets

typename = {'Cycles';'FWHM'};

figure(1), clf

for typei=1:2
    subplot(2,3,1+(typei-1)*3)
    contourf(timev,frex,real( squeeze(wavefam(typei,:,:)) ),40,'linecolor','none')
    set(gca,'clim',[-1 1]), axis square
    title([ typename{typei} ': real part' ])
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    
    subplot(2,3,2+(typei-1)*3)
    contourf(timev,frex,imag( squeeze(wavefam(typei,:,:)) ),40,'linecolor','none')
    set(gca,'clim',[-1 1]), axis square
    title([ typename{typei} ': imag part' ])
    
    subplot(2,3,3+(typei-1)*3)
    contourf(timev,frex,abs( squeeze(wavefam(typei,:,:)) ),40,'linecolor','none')
    set(gca,'clim',[-1 1]), axis square
    title([ typename{typei} ': magnitude' ])
end


%% show an example of one wavelet

% frex index
x = 99;

figure(2), clf
subplot(211), hold on
plot(timev,squeeze(real(wavefam(1,x,:))),'b')
plot(timev,squeeze(imag(wavefam(1,x,:))),'r')
plot(timev,squeeze(abs(wavefam(1,x,:))),'k')
legend({'real';'imag';'abs'})
title([ 'Cycles wavelet at ' num2str(frex(x)) ' Hz' ])

subplot(212), hold on
plot(timev,squeeze(real(wavefam(2,x,:))),'b')
plot(timev,squeeze(imag(wavefam(2,x,:))),'r')
plot(timev,squeeze(abs(wavefam(2,x,:))),'k')
legend({'real';'imag';'abs'})
title([ 'FWHM wavelet at ' num2str(frex(x)) ' Hz' ])
xlabel('Time (s)')

%% done.
