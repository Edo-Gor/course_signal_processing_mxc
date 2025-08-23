%%
%     COURSE: Solved problems in neural time series analysis
%    SECTION: Time-domain analyses
%      VIDEO: Butterfly plot and topo-variance time series
% Instructor: sincxpress.com
%
%%

load sampleEEGdata.mat

figure(1), clf

% make a butterfly plot (i.e., all ERP from all 64 channels overlayed)
subplot(211)
plot(EEG.times,squeeze(mean(EEG.data,3)),'linew',1.5)
set(gca,'xlim',[-500 1300])
title('Butterfly plot')
xlabel('Time (s)'), ylabel('Voltage (\muV)')
grid on
title('Earlobe reference')


%% compute the average reference
% this is the actual reason wy this is called 'butterfly' plot; being a
% reference to the mean, if summed together, the signal would be zero

% the fancy way...
data = bsxfun(@minus,EEG.data,mean(EEG.data,1));


subplot(212)
plot(EEG.times,squeeze(mean(data,3)),'linew',1.5)
set(gca,'xlim',[-500 1300])
title('Butterfly plot')
xlabel('Time (s)'), ylabel('Voltage (\muV)')
grid on
title('Average reference')
hold off

%% compute the variance time series

figure(2), clf

% variance for both earlobe and average reference; useful to have some more
% information about the above plots (indeed, as the mean is zero, what is 
% interesting is) the variance pattern)
var_ts_ER = var( mean(EEG.data,3) );  % earlobe refs; var over channels
var_ts_AR = var( mean(data,3) );      % average ref; var over channels

% of course the variance time series is the same as its computation
% involves basically an average referencing; tihs wouldn't be the same for
% the root mean square, which doesn't involve an average centering
subplot(211), hold on
plot(EEG.times,var_ts_ER,'s-','linew',1.5)
plot(EEG.times,var_ts_AR,'-','linew',1.5)

set(gca,'xlim',[-500 1300])
xlabel('Time (s)'), ylabel('Voltage (\muV)')
grid on
title('Topographical variance time series')

%% compute the root-mean-square

rms_ts_ER = rms( mean(EEG.data,3) );  % earlobe refs; rms over channels
rms_ts_AR = rms( mean(data,3) );      % average ref; rms over channels

subplot(212), hold on
plot(EEG.times,rms_ts_ER,'s-','linew',1.5)
plot(EEG.times,rms_ts_AR,'-','linew',1.5);

set(gca,'xlim',[-500 1300])
xlabel('Time (s)'), ylabel('Voltage (\muV)')
grid on
title('Topographical root-mean-square time series')
zoom on
hold off

%% done.
