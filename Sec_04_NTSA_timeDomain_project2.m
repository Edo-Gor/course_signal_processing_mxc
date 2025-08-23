%%
%     COURSE: Solved problems in neural time series analysis
%    SECTION: Time-domain analyses
%      VIDEO: Project 2-2: ERP peak latency topoplot
% Instructor: sincxpress.com
%
%%

%% Loop through each channel and find the peak time of the ERP between 100 and 400 ms. 
%   Store these peak times in a separate variable, and then make a
%   topographical plot of the peak times. Repeat for a low-pass filtered ERP.

load sampleEEGdata.mat

% time window 
boundaries = [100 400];
maxpeaktime = dsearchn(EEG.times',boundaries')';

% compute ERP and plot
erp = mean(EEG.data,3);

figure(1), clf
plot(EEG.times,erp)
set(gca,'xlim',[-500 1300])
title('Butterfly plot (earlobe ref. only)')
xlabel('Time (s)'), ylabel('Voltage (\muV)')
grid on

ylim = get(gca,'ylim');
ph = patch(EEG.times(maxpeaktime([1 2 2 1])),ylim([1 1 2 2]),'b');
set(ph,'facealpha',.2,'edgecolor','none')
set(gca,'Children',flipud(get(gca,'Children')))

% find max(erp) timing indices for each channel
[erpMax,erpMaxTime] = deal(zeros(1,EEG.nbchan));
for i = 1:EEG.nbchan
    [erpMax(i),erpMaxTime(i)] = max(erp(i,maxpeaktime(1):maxpeaktime(2)));
end

% bring indices into time
erpMaxTime_t = EEG.times(maxpeaktime(1) + erpMaxTime - 1);

% topoplot
figure(2), clf
subplot(121)
topoplotIndie(erpMaxTime_t, EEG.chanlocs,'numcontour',10,'electrodes','numbers','shading','interp');
set(gca,'clim',[100 400])
colorbar('eastoutside')
title({'ERP peak time';[' (' num2str(boundaries(1)) ' - ' num2str(boundaries(2)) ' ms)' ]})


%% Compute low-pass filter

% compute filter (15 Hz to reproduce example)
% set filter parameters
lowcut = 15;
filttime = -.3:1/EEG.srate:.3;
filtkern = sin(2*pi*lowcut*filttime) ./ filttime;

% adjust NaN and normalize filter to unit-gain
filtkern(~isfinite(filtkern)) = max(filtkern);
filtkern = filtkern./sum(filtkern);

% windowed sinc filter
filtkern = filtkern .* hann(length(filttime))';

% inspect the filter kernel
figure(3), clf
subplot(311)
plot(filttime,filtkern,'k','linew',1.5)
xlabel('Time (s)')
title('Time domain')

subplot(312)
hz = linspace(0,EEG.srate,length(filtkern));
plot(hz,abs(fft(filtkern)).^2,'ks-','linew',1.5)
set(gca,'xlim',[0 lowcut*3])
xlabel('Frequency (Hz)'), ylabel('Gain')
title('Frequency domain')

subplot(313)
hz = linspace(0,EEG.srate,length(filtkern));
plot(hz,abs(fft(filtkern)).^2,'ks-','linew',1.5)
set(gca,'xlim',[0 lowcut*3])
xlabel('Frequency (Hz)'), ylabel('Gain')
title('Frequency domain')
set(gca,'yscale','log')  % check also log scale

%% Low-pass filtered ERP

% compute filtered ERP
erp = double(erp);
erpF = zeros(size(erp));
for i = 1:EEG.nbchan
    erpF(i,:) = filtfilt(filtkern,1,erp(i,:));
end

% find max(erp) timing indices for each channel
[erpMaxF,erpMaxTimeF] = deal(zeros(1,EEG.nbchan));
for i = 1:EEG.nbchan
    [erpMaxF(i),erpMaxTimeF(i)] = max(erpF(i,maxpeaktime(1):maxpeaktime(2)));
end

% bring indices into time
erpMaxTimeF_t = EEG.times(maxpeaktime(1) + erpMaxTimeF - 1);

% topoplot
figure(2), hold on
subplot(122)
topoplotIndie(erpMaxTimeF_t, EEG.chanlocs,'numcontour',10,'electrodes','numbers','shading','interp');
set(gca,'clim',[100 400])
colorbar('eastoutside')
title({'Filtered ERP peak time';[' (' num2str(boundaries(1)) ' - ' num2str(boundaries(2)) ' ms)' ]})

%% Inspect differences 
% check e.g. channel 46, where the filtered peack is shifted because the
% unfiltered one gets considerably smoothed down 
chan_n = 47;

figure(4), clf
y_unf = erpMax(chan_n);
y_fil = erpMaxF(chan_n);
plot(EEG.times,erp(chan_n,:) ,EEG.times,erpF(chan_n,:))
yline(y_unf,'k--'); yline(y_fil,'k--')
set(gca,'xlim',boundaries)
xlabel(['Time window (' num2str(boundaries(1)) ' - ' num2str(boundaries(2)) ' ms)']), ylabel('Voltage (\muV)')
title(['A closer look to ERPs in channel ' num2str(chan_n) '.'])
zoom on



