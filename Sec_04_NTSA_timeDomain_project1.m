%%
%     COURSE: Solved problems in neural time series analysis
%    SECTION: Time-domain analyses
%      VIDEO: Project 2-1: Quantify the ERP as peak-mean or peak-to-peak
% Instructor: sincxpress.com
%
%%

load sampleEEGdata.mat

% channel to pick
chan2use = 'o1';
chandix = strcmpi({EEG.chanlocs.labels},chan2use);

% time window for negative peak
% Mike's choices: 
negpeaktime = [  50 110 ];
pospeaktime = [ 110 170 ];


%%% compute ERP
erp = double( mean(EEG.data(chandix,:,:),3) );

% plot ERP
figure(1), clf
plot(EEG.times,erp,'k','linew',1)
set(gca,'xlim',[-300 1000])

% plot patches over areas
ylim = get(gca,'ylim');
ph1 = patch(negpeaktime([1 2 2 1]),ylim([1 1 2 2]),'y');
ph2 = patch(pospeaktime([1 2 2 1]),ylim([1 1 2 2]),'g');
set(ph1,'FaceAlpha',.7,'EdgeColor','none')
set(ph2,'FaceAlpha',.7,'EdgeColor','none')


% move the patches to the background
set(gca,'Children',flipud( get(gca,'Children') ))

% titles 
xlabel('Times (ms)');
ylabel('Voltage (\muV)')
title(['ERP from channel ' chan2use])

%% first low-pass filter (windowed sinc function)

lowcut = 15;
filttime = -.3:1/EEG.srate:.3;
filtkern = sin(2*pi*lowcut*filttime) ./ filttime;

% adjust NaN and normalize filter to unit-gain
filtkern(~isfinite(filtkern)) = max(filtkern);
filtkern = filtkern./sum(filtkern);

% windowed sinc filter
filtkern = filtkern .* hann(length(filttime))';


% inspect the filter kernel
figure(2), clf
subplot(211)
plot(filttime,filtkern,'k','linew',2)
xlabel('Time (s)')
title('Time domain')


subplot(212)
hz = linspace(0,EEG.srate,length(filtkern));
spectr = abs(fft(filtkern)).^2;
plot(hz,spectr,'ks-','linew',2)
set(gca,'xlim',[0 lowcut*3])
xlabel('Frequency (Hz)'), ylabel('Gain')
title('Frequency domain')
% set(gca,'yscale','log')  % if you want log scale

%% now filter the ERP and replot

% apply filter (mind that it takes double as input)
erp_f = filtfilt(filtkern,1,erp);

% plot on top of unfiltered ERP
figure(1), hold on
plot(EEG.times,erp_f,'r','linew',1)
set(gca,'xlim',[-300 1000])
hold off

%% peak-to-peak voltages and timings

%%%% first for unfiltered ERP

% find index range of the 2 time windows
negpeak_i = dsearchn(EEG.times',negpeaktime');
pospeak_i = dsearchn(EEG.times',pospeaktime');

% find the min/max_peak in time-wind_1/2 (arg 1) and their inices WITHIN
% the respective time windows (arg 2)
[negpeak, neg_i] = min(erp(negpeak_i(1):negpeak_i(2)));
[pospeak, pos_i] = max(erp(pospeak_i(1):pospeak_i(2)));

% find index of min/max_peak within their respective windows
negpeak_t = negpeak_i(1) + neg_i - 1;
pospeak_t = pospeak_i(1) + pos_i - 1;

% compute peak-to-peak voltage and latency
erpP2P = abs(negpeak) + abs(pospeak);
erpP2Plat = EEG.times(pospeak_t) - EEG.times(negpeak_t);


%%%% then for low-pass filtered ERP

% find index range of the 2 time windows (not necessary to re-compute
% actually, just for clarity)
negpeak_i_F = dsearchn(EEG.times',negpeaktime');
pospeak_i_F = dsearchn(EEG.times',pospeaktime');

% find the min/max_peak in time-wind_1/2 (arg 1) and their inices WITHIN
% the respective time windows (arg 2)
[negpeak_F, neg_i_F] = min(erp_f(negpeak_i_F(1):negpeak_i_F(2)));
[pospeak_F, pos_i_F] = max(erp_f(pospeak_i_F(1):pospeak_i_F(2)));

% find index of min/max_peak within their respective windows
negpeak_t_F = negpeak_i_F(1) + neg_i_F -1;
pospeak_t_F = pospeak_i_F(1) + pos_i_F -1;

% compute peak-to-peak voltage and latency
erpFP2P = abs(negpeak_F) + abs(pospeak_F);
erpFP2Plat = EEG.times(pospeak_t_F) - EEG.times(negpeak_t_F);


%% Report the results in the command window

% clear the screen
clc

fprintf('\nRESULTS FOR PEAK POINT:')
fprintf('\n   Peak-to-peak on unfiltered ERP: %5.4g muV, %4.3g ms span.',erpP2P,erpP2Plat)
fprintf('\n   Peak-to-peak on filtered ERP:   %5.4g muV, %4.3g ms span.\n\n',erpFP2P,erpFP2Plat)

% results:
% RESULTS FOR PEAK POINT:
%    Peak-to-peak on unfiltered ERP: 11.27 muV, 66.4 ms span.
%    Peak-to-peak on filtered ERP:    8.05 muV, 58.6 ms span.

%%

%% repeat for mean around the peak

% time window for averaging (one-sided!!)
win = 10; % in ms

% now convert to indices (i.e. how many timepoints indeces can fit inside
% this thime window given our sampling rate)
winidx = round(win/(1000/EEG.srate));

%%%% first for unfiltered ERP

% compute the indices for the 2 times windows
negpeak_win2ave = negpeak_t-winidx : negpeak_t+winidx;
pospeak_win2ave = pospeak_t-winidx : pospeak_t+winidx;

% compute the average around the 2 peaks and then the mean(peak-to-peak)
% voltage; compute also the duration (same as above actually)
erpP2P = abs(mean(erp(negpeak_win2ave))) + abs(mean(erp(pospeak_win2ave)));


%%%% then for low-pass filtered ERP

% compute the indices for the 2 times windows
negpeak_win2ave_F = negpeak_t_F-winidx : negpeak_t_F+winidx;
pospeak_win2ave_F = pospeak_t_F-winidx : pospeak_t_F+winidx;

% compute the average around the 2 peaks and then the mean(peak-to-peak)
% voltage; compute also the duration (same as above actually)
erpFP2P = abs(mean(erp_f(negpeak_win2ave_F))) + abs(mean(erp_f(pospeak_win2ave_F)));


%% Report the results in the command window

fprintf('\nRESULTS FOR WINDOW AROUND PEAK:')
fprintf('\n   Peak-to-peak on unfiltered ERP: %5.4g muV, %4.3g ms span.',erpP2P,erpP2Plat)
fprintf('\n   Peak-to-peak on filtered ERP:   %5.4g muV, %4.3g ms span.\n\n',erpFP2P,erpFP2Plat)

% results:
% RESULTS FOR WINDOW AROUND PEAK:
%    Peak-to-peak on unfiltered ERP: 7.996 muV, 66.4 ms span.
%    Peak-to-peak on filtered ERP:   7.442 muV, 58.6 ms span.


%% done.
