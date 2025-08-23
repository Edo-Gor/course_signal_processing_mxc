%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Connectivity analyses
%      VIDEO: Power correlations over trials
% Instructor: sincxpress.com
%
%%

clear

% load EEG data
load sampleEEGdata.mat

% time-frequency analysis parameters
freqrange = [2 25]; % frequency range in Hz
numfrex   = 30;
frex      = linspace(freqrange(1),freqrange(2),numfrex);
fwhms     = linspace(.5,.2,numfrex);
baseidx   = dsearchn(EEG.times',-[500 200]');


% connectivity analysis parameters
tidx = dsearchn(EEG.times',[-600 -200]');
fidx = dsearchn(frex',[14 20]');

% channel to use
channel = strcmpi('fcz',{EEG.chanlocs.labels});

%% setup time-frequency analysis (parameters and wavelets)

% other convolution parameters
frex    = linspace(freqrange(1),freqrange(2),numfrex);
wavtime = -2:1/EEG.srate:2;
nData   = EEG.pnts*EEG.trials;
nKern   = length(wavtime);
nConv   = nData + nKern - 1;
halfwav = (length(wavtime)-1)/2;

% create wavelets
cmwX = zeros(numfrex,nConv);
for fi=1:numfrex
    
    % create time-domain wavelet
    gausWin = exp( -4*log(2)*wavtime.^2 / fwhms(fi)^2 );
    cmw = exp(2*1i*pi*frex(fi).*wavtime) .* gausWin;
    
    % compute fourier coefficients of wavelet and normalize
    cmwX(fi,:) = fft(cmw,nConv);
    cmwX(fi,:) = cmwX(fi,:) ./ max(cmwX(fi,:));
end

%% time-frequency analysis

% initialize time-frequency output matrices
tf = zeros(numfrex,EEG.pnts);
tf3d = zeros(numfrex,EEG.pnts,EEG.trials);

% compute Fourier coefficients of EEG data (doesn't change over frequency!)
eegX = fft( reshape(EEG.data(channel,:,:),1,[]) ,nConv);

% loop over frequencies
for fi=1:numfrex
    
    % second and third steps of convolution
    as = ifft( eegX .* cmwX(fi,:) );
    
    % cut wavelet back to size of data
    as = as(halfwav+1:end-halfwav);
    as = reshape(as,EEG.pnts,EEG.trials);
    
    % trial-averaged power
    powerts = mean(abs(as).^2,2);
    tf(fi,:) = 10*log10( powerts/mean(powerts(baseidx(1):baseidx(2))) );
    
    % extract power from all trials
    tf3d(fi,:,:) = abs(as).^2;
    
end % end frequency loop

% reshape the 3D matrix to 2D
tf2d = reshape(tf3d,numfrex*EEG.pnts,EEG.trials)';

%% show TF power

figure(1), clf
contourf(EEG.times,frex,tf,40,'linecolor','none')
set(gca,'clim',[-3 3],'xlim',[-300 1200])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Trial-averaged time-frequency power from ' EEG.chanlocs(channel).labels ])
colormap jet

%% extract seed power and create design matrix

% extract power on each trial
seed_power = squeeze(mean(mean(tf3d(fidx(1):fidx(2),tidx(1):tidx(2),:),1),2));

% create new design matrix with power
designMat = [ ones(EEG.trials,1) seed_power ];

%% show the design and data matrices

figure(2), clf

ax1_h = axes;
set(ax1_h,'Position',[.05 .1 .1 .8])
imagesc(designMat)
set(ax1_h,'xtick',1:2,'xticklabel',{'Int';'RTs'},'ydir','norm')
ylabel('Trials')
title('Design matrix')


ax2_h = axes;
set(ax2_h,'Position',[.25 .1 .7 .8])
imagesc(tf2d)
set(ax2_h,'ydir','norm','clim',[0 20])
ylabel('Trials')
xlabel('Timefrequency')
title('Data matrix')

colormap gray

%% connectivity using least-squares modeling

% least-squares model
x = (designMat'*designMat)\designMat'*tf2d;

% extract only the second regressor (seed power)
betamat = reshape(x(2,:),numfrex,EEG.pnts);


% now show the results
figure(3), clf
% show time-frequency map of regressors
contourf(EEG.times,frex,betamat,40,'linecolor','none')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'xlim',[-700 1200],'clim',[-2 2])

% box to show 'seed' region
hold on
plot(EEG.times(tidx),[frex(fidx(1)) frex(fidx(1))],'k','linew',2)
plot(EEG.times(tidx),[frex(fidx(2)) frex(fidx(2))],'k','linew',2)
plot([EEG.times(tidx(1)) EEG.times(tidx(1))],frex(fidx),'k','linew',2)
plot([EEG.times(tidx(2)) EEG.times(tidx(2))],frex(fidx),'k','linew',2)

title([ 'Correlations with power from ' num2str(round(EEG.times(tidx(1)))) ' to ' num2str(round(EEG.times(tidx(2)))) ', ' num2str(round(frex(fidx(1)))) '-' num2str(round(frex(fidx(2)))) ' Hz' ])
colormap jet

%% repeat with simple correlations

% initialize correlation matrix
cormat = zeros(numfrex,EEG.pnts);

% loop over all time-frequency points
for fi=1:numfrex
    for ti=1:EEG.pnts
        
        % correlation at this TF point
        tfpoint = tf3d(fi,ti,:);
        cormat(fi,ti) = corr(seed_power,tfpoint(:));
        
    end
end


figure(4), clf
% show time-frequency map of regressors
contourf(EEG.times,frex,cormat,40,'linecolor','none')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'xlim',[-700 1200],'clim',[-1 1]*.7)

% box to show 'seed' region
hold on
plot(EEG.times(tidx),[frex(fidx(1)) frex(fidx(1))],'k','linew',2)
plot(EEG.times(tidx),[frex(fidx(2)) frex(fidx(2))],'k','linew',2)
plot([EEG.times(tidx(1)) EEG.times(tidx(1))],frex(fidx),'k','linew',2)
plot([EEG.times(tidx(2)) EEG.times(tidx(2))],frex(fidx),'k','linew',2)

title([ 'Correlations with power from ' num2str(round(EEG.times(tidx(1)))) ' to ' num2str(round(EEG.times(tidx(2)))) ', ' num2str(round(frex(fidx(1)))) '-' num2str(round(frex(fidx(2)))) ' Hz' ])
colormap jet 

%% done.
