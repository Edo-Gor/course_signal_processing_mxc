%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Connectivity
%      VIDEO: Spurious connectivity in narrowband noise
% Instructor: sincxpress.com
%
%%

clear

% simulation details
pnts  = 10000;
srate =  1000;
time  = (0:pnts-1)/srate;


% create filtered noise
frange = [10 15]; % in frequency
order  = round( 20*srate/frange(1) );

% filter kernel
filtkern = fir1(order,frange/(srate/2));

% signals are filtered noise
signal1 = filtfilt(filtkern,1,randn(1,pnts));
signal2 = filtfilt(filtkern,1,randn(1,pnts));


% plot the signals
figure(1), clf
subplot(311)
plot(time,signal1, time,signal2)
ylabel('Amplitude')
title('Two signals in time domain')

%% calculate phase synchronization

% extract angles from Hilbert transform
angles1 = angle(hilbert( signal1 ));
angles2 = angle(hilbert( signal2 ));

% show phase angle time series
subplot(312)
plot(time,angles1, time,angles2)
ylabel('Phase angle (rad.)')
title('Phase angle time series')

%% synchronization over time (in windows)

% analysis parameters
winsize = 300; % in points

% initialize
synchTS = zeros(1,pnts);

% loop over time points
for ti=1:pnts
    
    % figure out time points, considering edges
    tidx = max(1,ti-winsize) : min(pnts,ti+winsize);
    
    % compute synchronization for this center time points
    synchTS(ti) = abs(mean(exp(1i*( angles1(tidx)-angles2(tidx) ))));
end


% and plot
subplot(313)
plot(time,synchTS,'k','linew',2)
set(gca,'ylim',[0 1])
xlabel('Time (s)'), ylabel('Synch. strength')
title('Synchronization over time')

%% done.
