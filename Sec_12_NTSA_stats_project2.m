%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Inferrential statistics
%      VIDEO: Project 7-2: Simulate time-frequency data for statistical testing
% Instructor: sincxpress.com
%
%%

clear

%% instructions:

%%% data:
% 2 signals, each with two non-phase-locked time-frequency features
% some differences between conditions
% 100 trials / condition

%%% noise:
% 1/f, trial-unique

%%% data analysis:
% wavelet convolution

%%% statistics
% mapwise permutation testing
% cluster correction
% separately for positive vs. negative features!


%% generic simulation parameters

npnts   = 2501;
time    = linspace(-1,3,npnts)'; % force column vector
srate   = 1./mean(diff(time));   % Hz (mean of diff(time) to account for instabilities)

ntrials = 100;                   % per condition
nPerms  = 500;                   % number of iterations in permutation testing
pval    = .05;                   % p-value threshold

sigThresh = norminv(1-pval/2);   % two-tailed!


% dataset parameters
freq1a = 6;
freq1b = 15;
freq2a = 7;
freq2b = 14;

peak1a = 1;
peak1b = 1.6;
peak2a = .9;
peak2b = 1.7;

amp1a  = 4;
amp1b  = 2;
amp2a  = 2;
amp2b  = 5;

fwhm1a = .5;
fwhm1b = .3;
fwhm2a = .4;
fwhm2b = .5;

noisestd = 1.5;

% create Gaussians for signal envelopes
gausw1a = exp( -4*log(2)*(time-peak1a).^2 / fwhm1a.^2 );
gausw1b = exp( -4*log(2)*(time-peak1b).^2 / fwhm1b.^2 );
gausw2a = exp( -4*log(2)*(time-peak2a).^2 / fwhm2a.^2 );
gausw2b = exp( -4*log(2)*(time-peak2b).^2 / fwhm2b.^2 );

% 1/f parameters
ed = 50; % exponential decay parameter

%% create the datasets
% loop to assure some new random noise at each trial

[data1,data2] = deal( zeros(npnts,ntrials) ); % time by trial
for triali=1:ntrials
    
    %%% dataset 1
    % create signal1
    sig1a = amp1a * cos(2*pi*freq1a*time + rand*2*pi) .* gausw1a;
    sig1b = amp1b * cos(2*pi*freq1b*time + rand*2*pi) .* gausw1b;
    
    % create 1/f noise (check scripts from section 3; basically create
    % spectrum and then compute ifft)
    as = rand(1,floor(npnts/2)-1) .* exp(-(1:floor(npnts/2)-1)/ed);
    as = [as(1) as 0 0 as(:,end:-1:1)]';
    fc = as .* exp(1i*2*pi*rand(size(as)));
    noise = real(ifft(fc));
    noise = noise*(noisestd/std(noise)); % adjust STD (check previous project)
    
    data1(:,triali) = sig1a + sig1b + noise;
    
    
    %%% dataset 2
    % create signal2
    sig2a = amp2a * cos(2*pi*freq2a*time + rand*2*pi) .* gausw2a;
    sig2b = amp2b * cos(2*pi*freq2b*time + rand*2*pi) .* gausw2b;
    
    % create 1/f noise
    as = rand(1,floor(npnts/2)-1) .* exp(-(1:floor(npnts/2)-1)/ed);
    as = [as(1) as 0 0 as(:,end:-1:1)]';
    fc = as .* exp(1i*2*pi*rand(size(as)));
    noise = real(ifft(fc));
    noise = noise*(noisestd/std(noise)); % adjust STD (check previous project)
    
    data2(:,triali) = sig2a + sig2b + noise;
    
end

% concatenate into one big matrix and generate condition vector for later 
dataFull = cat(2,data1,data2);
condlabels = (1:ntrials*2)>ntrials;


% plot
figure(1), clf
plot(time,mean(data1,2), time,mean(data2,2))
ylabel('Amplitude (a.u.)')
legend({'condition 1';'condition 2'})
title([ 'Signals plus noise (average of ' num2str(ntrials) ' trials)' ])


%% time-frequency analysis via wavelet convolution

% frequency parameters
minFreq =  2;
maxFreq = 30;
numFrex = 50;
frex = linspace(minFreq,maxFreq,numFrex);

wavefwhms = linspace(1/minFreq,2/maxFreq,numFrex);
wavtime = -2:1/srate:2;
half_wave = (length(wavtime)-1)/2;

% FFT parameters
nWave = length(wavtime);
nData = npnts * 2*ntrials;
nConv = nWave + nData - 1;


% now compute the FFT of all trials concatenated
dataX = fft( reshape( dataFull,1,[]) ,nConv );


% initialize output time-frequency data and perform convolution, looping 
% over frequencies
tf = zeros(numFrex,npnts,2*ntrials); 
for fi=1:length(frex)
    
    % create wavelet and get its FFT
    wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*wavefwhms(fi)^2));
    waveletX = fft(wavelet,nConv);
    waveletX = waveletX ./ max(waveletX); % max normalise
    
    % now run convolution in one step (as product in freq. domain)
    as = ifft(waveletX .* dataX);
    as = as(half_wave+1:end-half_wave);
    
    % and reshape back to time X trials
    as = reshape( as, npnts, 2*ntrials );
    
    % compute power and average over trials
    tf(fi,:,:) = abs(as).^2;
end


% plot the time-frequency results
figure(2), clf
for i=1:2
    subplot(2,1,i)
    imagesc(time,frex,squeeze(mean(tf(:,:,condlabels==i-1),3)))
    axis xy
    set(gca,'clim',[0 .1])
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    title([ 'Condition ' num2str(i) ])
    colormap jet 
end


%% permutation testing

% initialize and compute pixel-level permutations (N=500)
permutedDiffs = zeros(nPerms,numFrex,npnts);
for permi = 1:nPerms
    
    % shuffle condition label vector
    fakeconds = condlabels( randperm(2*ntrials) );
    
    % compute and store difference time series (average over the time
    % domain)
    mean1 = mean( tf(:,:,fakeconds==0),3 );
    mean2 = mean( tf(:,:,fakeconds==1),3 );
    permutedDiffs(permi,:,:) = mean1-mean2;
end


% examine differences without cluster correction with some plotting before 
% cluster correction is applied

% difference between 'true' data (cond. 1 and 2; time by frequencies)
obsdiff = mean(tf(:,:,condlabels==0),3) - mean(tf(:,:,condlabels==1),3);
% compute z maps (time by frequencies); that is, the difference between the
% true data and the average of the permutations
zmap = (obsdiff - squeeze(mean(permutedDiffs,1))) ./ squeeze(std(permutedDiffs));

figure(3), clf
% raw difference map
subplot(221)
contourf(time,frex,obsdiff,40,'linecolor','none')
set(gca,'clim',[-1 1]*.5)
title('Difference map')

% difference z-map
subplot(222)
contourf(time,frex,zmap,40,'linecolor','none')
set(gca,'clim',[-1 1]*7)
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Difference z-map')
colormap jet 

% thresholded zmap (zero out values less than sigThresh; still pixel-wise)
zthresh = zmap;
zthresh( abs(zmap)<sigThresh ) = 0;

subplot(223)
contourf(time,frex,zthresh,40,'linecolor','none')
set(gca,'clim',[-1 1]*7)
title('Thresholded difference z-map')

% raw difference map with contours
subplot(224), hold on
contourf(time,frex,obsdiff,40,'linecolor','none')
contour(time,frex,logical(zthresh),1,'linecolor','k')
set(gca,'clim',[-1 1]*.5)
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Difference map with thresholded outlines')

%% find cluster sizes under the null hypothesis

% initialize cluster sizes from permutation
clustsizes = zeros(nPerms,2);

% precompute mean and std images
permMean = squeeze( mean(permutedDiffs) );
permStd  = squeeze( std(permutedDiffs) );

for permi=1:nPerms
    
    % POSITIVE-valued clusters
    zdiffFake = (squeeze(permutedDiffs(permi,:,:))-permMean) ./ permStd;
    zdiffFake( zdiffFake<sigThresh ) = 0;
    islands = bwconncomp( logical(zdiffFake) );
    clustNs = cellfun(@length,islands.PixelIdxList);
    clustsizes(permi,1) = max(clustNs);
    
    % NEGATIVE-valued clusters
    zdiffFake = (squeeze(permutedDiffs(permi,:,:))-permMean) ./ permStd;
    zdiffFake( zdiffFake>-sigThresh ) = 0;
    islands = bwconncomp( logical(zdiffFake) );
    clustNs = cellfun(@length,islands.PixelIdxList);
    clustsizes(permi,2) = max(clustNs);
end

% compute EMPIRICAL cluster threshold
clustthresh(1) = prctile(clustsizes(:,1),100-(pval/2)*100); % positive clusters 
clustthresh(2) = prctile(clustsizes(:,2),100-(pval/2)*100); % negative clusters


% show distribution of cluster sizes
figure(4), clf
for i=1:2
    subplot(2,1,i)
    histogram(clustsizes(:,i),30)
    xlabel('Cluster size (pixels)'), ylabel('Count')
    if i==1
        title('Positive clusters') 
    else
        title('Negative clusters')
    end
end

%% remove small clusters from real thresholded data

% thresholded zmap (separate for each tail)
[zthreshP,zthreshN] = deal( zmap );
zthreshP( zmap<sigThresh ) = 0;     % zero out positive tail
zthreshN( zmap>-sigThresh ) = 0;    % zero out negative tail

% find islands (separately)
islandsP = bwconncomp( logical(zthreshP) );
islandsN = bwconncomp( logical(zthreshN) );

% find and remove any subthreshold islands
for ii=1:islandsP.NumObjects
    if numel(islandsP.PixelIdxList{ii})<clustthresh(1)
        zthreshP(islandsP.PixelIdxList{ii}) = 0;
    end
end
% repeat for negative
for ii=1:islandsN.NumObjects
    if numel(islandsN.PixelIdxList{ii})<clustthresh(2)
        zthreshN(islandsN.PixelIdxList{ii}) = 0;
    end
end


% now plot it
figure(5), clf, hold on
contourf(time,frex,zmap,40,'linecolor','none')
contour(time,frex,logical(zthreshP),1,'linecolor','k','linew',2)
contour(time,frex,logical(zthreshN),1,'linecolor','r','linew',2)
set(gca,'clim',[-1 1]*9)
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Difference map with thresholded outlines, cluster-corrected')
colorbar
colormap jet 






