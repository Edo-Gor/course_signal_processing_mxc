%%
%     COURSE: Solved problems in neural time series analysis
%    SECTION: Time-frequency analyses
%      VIDEO: Project 4-1: Phase-locked, non-phase-locked, and total power
% Instructor: sincxpress.com
%
%%

clear, clc

%% 

load v1_laminar.mat

chanidx = 7;

% wavelet parameters
minFreq = 1;
maxFreq = 80;
nFrex   = 76;
frex    = linspace(minFreq,maxFreq,nFrex);

% other wavelet parameters
wtime   = -1:1/srate:1-1/srate;
halfwav = floor(length(wtime)/2);
FWHMrange = [.4 .2];
FWHMs     = logspace(log10(FWHMrange(1)),log10(FWHMrange(2)),nFrex);

% baseline time window
baseline_time = [ -.4 -.1 ];

% FFT parameters
n_wave = length(wtime);
n_data = length(timevec)*size(csd,3);
n_conv = n_wave + n_data -1;

%% create non-phase-locked dataset

% compute ERP
erp = mean(csd(chanidx,:,:),3);

% compute induced power by subtracting ERP from each trial
nonphaselocked = squeeze(csd(chanidx,:,:)) - repmat(erp',1,size(csd,3));

% FFT of data
dataX{1} = fft(reshape(csd(chanidx,:,:),1,[]),n_conv);
dataX{2} = fft(reshape(nonphaselocked(:,:),1,[]),n_conv);

% convert baseline from ms to indices
bidx = dsearchn(timevec',baseline_time');

%% run convolution

% initialize output time-frequency data
tf = zeros(2,length(frex),length(timevec));

% loop over frequencies
for fi=1:length(frex)
    
    % create wavelet
    cmw = exp( 1i*2*pi*frex(fi)*wtime ) .* exp( -4*log(2)*wtime.^2/FWHMs(fi)^2 );
    % take FFT of wavelet and max scale (even though max scale not necessary)
    cmwX = fft(cmw,n_conv);
    cmwX = cmwX ./ max(cmwX);
    
    % run convolution for total and non-phase-locked
    for i=1:2
        
        % convolve, bring back to time domain and clip wings
        convres = ifft( cmwX.*dataX{i} );
        convres = convres(halfwav+1:end-halfwav);
        
        % reshape back to timeXtrials
        convres = reshape(convres,size(csd,2),size(csd,3));
        
        % compute power, averaged across trials
        tf(i,fi,:) = mean(abs(convres).^2,2);
        
        % db correct power
        basepow = mean(mean(abs(convres(bidx(1):bidx(2),:)).^2,1),2); 
        tf(i,fi,:) = 10*log10( tf(i,fi,:) / basepow );
        % or ...
%         tf(i,fi,:) = 10*log10( squeeze(tf(i,fi,:)) ./ mean(tf(i,fi,bidx(1):bidx(2)),3) );
        
    end % end loop around total/nonphase-locked/phaselocked
end % end frequency loop

%% plotting

analysis_labels = {'Total';'Non-phase-locked'};

% color limits
clim = [-15 15];

% scale ERP for plotting
erpt = (erp-min(erp))./max(erp-min(erp));
erpt = erpt*(frex(end)-frex(1))+frex(1);
erpt = erpt/2+10;

figure(1), clf
for i=1:2
    
    subplot(1,3,i), hold on
    contourf(timevec,frex,squeeze(tf(i,:,:)),40,'linecolor','none')
    
    set(gca,'clim',clim,'xlim',[-.2 1.5])
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
    title(analysis_labels{i})
    axis square
    colormap jet
    colorbar

    % plot ERP on top
    if i == 1
        plot(timevec,erpt,'k')
    end
end

% phase-locked component is the difference between total and non-phase-locked
subplot(133), hold on
contourf(timevec,frex,squeeze(tf(1,:,:)-tf(2,:,:)),40,'linecolor','none')
set(gca,'clim',clim/1,'xlim',[-.2 1.5])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Phase-locked')
colormap jet
colorbar
plot(timevec,erpt,'k')
axis square

%% end.
