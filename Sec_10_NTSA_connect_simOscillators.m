%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Connectivity
%      VIDEO: Synchronization in simulated noisy oscillators
% Instructor: sincxpress.com
%
%%

% simulation details
pnts  = 7000;
srate = 1000;
time  = (0:pnts-1)/srate;


% create filtered noise
frange = [10 15]; % in frequency
order  = round( 20*srate/frange(1) );

% filter kernel (narrowband filter)
filtkern = fir1(order,frange/(srate/2));

% create noise (filtered noise)
noise1 = 10 * filtfilt(filtkern,1,randn(1,pnts));
noise2 = 10 * filtfilt(filtkern,1,randn(1,pnts));

% signal is pure sine wave plus noise
signal1 = sin(2*pi*mean(frange)*time) + noise1;
signal2 = sin(2*pi*mean(frange)*time) + noise2;


% plot the signals
figure(1), clf
subplot(211)
plot(time,signal1, time,signal2)
xlabel('Time (s)'), ylabel('Amplitude')
title('Two signals in time domain')

%% calculate phase synchronization

% extract angles from Hilbert transform
angles1 = angle(hilbert( signal1 ));
angles2 = angle(hilbert( signal2 ));

% show phase angle time series
subplot(212)
plot(time,angles1, time,angles2)
xlabel('Time (s)'), ylabel('Phase angle (rad.)')
title('Phase angle time series')

% synchronization as average vector
avevec = mean(exp(1i*( angles1-angles2 )));

%% more plotting

figure(2), clf

% polar plot of angles from the two oscillators
subplot(121)
N2plot = 200;
h1 = polar([zeros(1,N2plot); angles1(1:N2plot)],[zeros(1,N2plot); ones(1,N2plot)]);
hold on
h2 = polar([zeros(1,N2plot); angles2(1:N2plot)],[zeros(1,N2plot); ones(1,N2plot)]);
% adjust the colors...
set(h1,'color','b')
set(h2,'color',ones(1,3)*.2)



% phase angle differences distribution
subplot(122)
% N2plot = pnts;
h = polar([zeros(1,N2plot); angles1(1:N2plot)-angles2(1:N2plot)],[zeros(1,N2plot); ones(1,N2plot)]);
set(h,'color',ones(1,3)*.5)
hold on

% now show the average phase angle...
h = polar([0 angle(avevec)],[0 abs(avevec)]);
set(h,'color','k','linew',3)

% ...and report the synchronization magnitude
title([ 'Phase synchronization: ' num2str( abs(avevec) ) ])

%% done.
