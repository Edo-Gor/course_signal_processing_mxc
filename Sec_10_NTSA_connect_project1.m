%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Connectivity
%      VIDEO: Project 5-1: ISPC and PLI, with and without Laplacian
% Instructor: sincxpress.com
%
%%

% load data and pick channels
clear
load sampleEEGdata

% compute Laplacian
EEG.lap = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);

%% select parameters

% pick two channels for synchronization
chan1 = 'FCz';
chan2 = 'POz';

% frequency range for convolution
min_freq =  2;
max_freq = 20;
num_frex = 30;
srate    = EEG.srate;

% set range for wavelet FWHM
% Note: encoded as number of cycles per frequency
fwhm_range = [ 2/min_freq 4/max_freq ];

% frequencies vectors
frex  = linspace(min_freq,max_freq,num_frex);
fwhms = linspace(fwhm_range(1),fwhm_range(2),num_frex);

%% wavelet parameters 

% wavelet and convolution parameters
wtime = -2:1/srate:2;
nWave = length(wtime);
nData = size(EEG.data,2)*size(EEG.data,3); % timepoints by trials
nConv = nData + nWave - 1;
halfw = (length(wtime)-1)/2;

%% FFT of data and Laplacian

% voltage data, chan1 and 2
EEG1X = fft( reshape(EEG.data(strcmpi(chan1,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv);
EEG2X = fft( reshape(EEG.data(strcmpi(chan2,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv );

% laplacian data, chan1 and 2
LAP1X = fft( reshape(EEG.lap(strcmpi(chan1,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv );
LAP2X = fft( reshape(EEG.lap(strcmpi(chan2,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv );


%% convolution etc

% initialize output time-frequency data
[ispc,pli] = deal( zeros(2,num_frex,EEG.pnts) );

% loop over frequencies
for fi=1:num_frex
    
    %%% create wavelet
    % create wavelet and get its FFT
    cmw  = exp( 1i*2*pi*frex(fi).*wtime ) .* exp( -4*log(2)*wtime.^2./fwhms(fi)^2 );
    cmwX = fft(cmw,nConv);
    % need to normalize?? nope, as far as we're looking for phases
    
    
    %%% convolution for voltage data
    % chan1
    as1 = ifft( EEG1X.*cmwX,nConv );
    as1 = reshape(as1(halfw+1:end-halfw),EEG.pnts,EEG.trials);
    
    % chan2
    as2 = ifft( EEG2X.*cmwX,nConv );
    as2 = reshape(as2(halfw+1:end-halfw),EEG.pnts,EEG.trials);
    
    % collect "eulerized" phase angle differences
    phasediffVOLT = exp(1i*( angle(as1)-angle(as2) ));
    
    
    %%% convolution for Laplacian data
    % chan1
    as1 = ifft( LAP1X.*cmwX,nConv );
    as1 = reshape(as1(halfw+1:end-halfw),EEG.pnts,EEG.trials);
    
    % chan2
    as2 = ifft( LAP2X.*cmwX,nConv );
    as2 = reshape(as2(halfw+1:end-halfw),EEG.pnts,EEG.trials);
    
    % collect "eulerized" phase angle differences
    phasediffLAP = exp(1i*( angle(as1)-angle(as2) ));
    
    
    %%% connectivities
    % ISPC and PLI for voltage
    ispc(1,fi,:) = abs(mean(phasediffVOLT,2));
    pli(1,fi,:)  = abs(mean(sign(imag(phasediffVOLT)),2));
    
    % ISPC and PLI for laplacian
    ispc(2,fi,:) = abs(mean(phasediffLAP,2));
    pli(2,fi,:)  = abs(mean(sign(imag(phasediffLAP)),2));
    
end


%% plotting

figure(1), clf
colormap jet

clim = [.1 .4];
datalabels = {'Voltage';'Laplacian'};

for i=1:2
    
    % ISPC
    subplot(2,2,i)
    contourf(EEG.times,frex,squeeze(ispc(i,:,:)),40,'linecolor','none')
    set(gca,'xlim',[-300 1000],'clim',clim)
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title([ 'ISPC: ' datalabels{i} ])
    
    
    % PLI
    subplot(2,2,i+2)
    contourf(EEG.times,frex,squeeze(pli(i,:,:)),40,'linecolor','none')
    set(gca,'xlim',[-300 1000],'clim',clim)
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title([ 'PLI: ' datalabels{i} ])
end

%% end.
