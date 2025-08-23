%%
%     COURSE: Solved problems in neural time series analysis
%    SECTION: Time-domain analyses
%      VIDEO: Lowpass filter an ERP
% Instructor: sincxpress.com
%
%%

load v1_laminar

% reduce data for convenience (channel 7)
data = double(squeeze( csd(7,:,:) ));


% cutoff frequency for low-pass filter
lowcut = 20; % in Hz

%% create and inspect the filter kernel

% filter order
order = 18; % notice how higher order filters have a sharper freq domain, but are also longer
filtord = round( order * (lowcut*1000/srate) );

% create filter
filtkern = fir1(filtord,lowcut/(srate/2),'low');

% inspect the filter kernel
figure(1), clf
subplot(211)
plot((0:length(filtkern)-1)/srate,filtkern,'k','linew',2)
xlabel('Time (s)')
title('Time domain')


subplot(212)
hz = linspace(0,srate,length(filtkern));
plot(hz,abs(fft(filtkern)).^2,'ks-','linew',2)
set(gca,'xlim',[0 lowcut*3])
xlabel('Frequency (Hz)'), ylabel('Gain')
title('Frequency domain')
set(gca,'yscale','log')

%% option 1: filter the ERP

% extract ERP
erp1 = mean(data,2); % dim 2 is trial here

figure(2) % quick plot to compare the un-filtered and the filtered version of erp
plot(timevec,erp1)
xlabel('Time (s)'), ylabel('Activity (\muV)')
axis([timevec(1),timevec(end),min(erp1)-50,max(erp1)+50])
hold on

% apply filter
erp1 = filtfilt(filtkern,1,erp1);
plot(timevec,erp1) % see how the time-locked high-freq activity at 0.7s is removed
hold off

%% option 2: filter the single trials instead of the overall ERP

erp2 = zeros(size(timevec));

for triali=1:size(data,2)
    erp2 = erp2 + filtfilt(filtkern,1,data(:,triali))';
end

% complete the averaging
erp2 = erp2/triali;

%% option 3: concatenate the trials, filter, and then separate them again

% make one long trial
supertrial = reshape(data,1,[]);

figure(3) % quick plot to compare the un-filtered and the filtered version of supertrial
timevec_scaled = linspace(timevec(1),range(timevec)*200,length(supertrial));
plot(timevec_scaled,supertrial)
xlabel('Time (s)'), ylabel('Activity (\muV)')
axis([timevec_scaled(1),timevec_scaled(end),min(supertrial)-50,max(supertrial)+50])
hold on

% apply filter
supertrial = filtfilt(filtkern,1,supertrial);
plot(timevec_scaled,supertrial) % zoom in to have a better look
hold off

% reshape back and take average
erp3 = reshape(supertrial,size(data));
erp3 = mean(erp3,2);

%% now compare

figure(4), clf, hold on

s = 'so^';

% cool way to loop over variables with different names (erp1, erp2, erp3
% ...), by concatenatig the end, instead of using plot 3 times
for i=1:3
    eval([ 'plot(timevec,erp' num2str(i) ',[s(i) ''-'' ],''linew'',1)' ])
end

axis([timevec(1),timevec(end),min(erp1)-50,max(erp1)+50])
xlabel('Time (s)'), ylabel('Voltage (\muV)')
legend({'filter ERP';'filter trials';'filter concat'})

% let's take a closer look
zoom on

% basically identical except for some small edge effect; it's a linear
% operation so it doesn't matter the order (e.g., filtering then averaging
% or viceversa); if we were to use non-linear operations, the order would
% have mattered

%% done.
