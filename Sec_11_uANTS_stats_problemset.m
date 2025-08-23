%%
%   COURSE: Neural signal processing and analysis: Zero to hero
%  SECTION: Problem set: Permutation-based statistics
%  TEACHER: Mike X Cohen, sincxpress.com
%


%% Did our little mousey get tired of this task? let's find out!
% Compute TF power, and store the single-trial power from the first 50 trials and last 50 trials.
% Then perform a statistical test of power (permutation testing) for early vs. late trials. 
% Use a p-value of .05. Show the results at an uncorrected level, and correcting for 
% multiple comparisons using cluster correction (or extreme-pixel if you don't have the image-processing toolbox). 
% Pick one electrode for this test.
% If you wanted to test the effects in three electrodes (e.g., deep, L4, superficial), 
%     what would be an appropriate p-value?


% This space intentionally left blank.

clear
load v1_laminar.mat

chan2use = 5;
times2save = -.3:.01:1.1;

xlim       = [-.1 1]; % for plotting

% specify frequency range
min_freq = 10;
max_freq = 100;
num_frex = 60;

% define frequency and time ranges
frex = linspace(min_freq,max_freq,num_frex);
times2saveidx = dsearchn(timevec',times2save');

% parameters for complex Morlet wavelets
wavtime  = -2:1/srate:2-1/srate;
half_wav = (length(wavtime)-1)/2;
cycRange = [ 4 12 ];
nCycles  = logspace(log10(cycRange(1)),log10(cycRange(end)),num_frex); % log10() for appropriate spacing

% FFT parameters
nWave = length(wavtime);
nData = size(csd,2)*size(csd,3);
nConv = nWave+nData-1;

% and create wavelets
cmwX = zeros(num_frex,nConv);
for fi=1:num_frex
    s          = nCycles(fi) / (2*pi*frex(fi));
    cmw        = exp(1i*2*pi*frex(fi).*wavtime) .* exp( (-wavtime.^2) ./ (2*s^2) );
    tempX      = fft(cmw,nConv);
    cmwX(fi,:) = tempX ./ max(tempX);
end

%% run convolution to extract tf power

% spectra of data; concatenate time and trials from channel of interest
dataX = fft( reshape(csd(chan2use,:,:),1,[]) ,nConv);

% initialize output time-frequency data
% notice that here we save all trials
tf = zeros(num_frex,length(times2save),100); % 100 because interested in the first 50 and 

for fi=1:num_frex
    
    % run convolution
    as = ifft(cmwX(fi,:).*dataX);
    as = as(half_wav+1:end-half_wav);           % cut wings
    as = reshape(as,size(csd,2),size(csd,3));   % resize
    
    % power on all trials from channel "1"
    % only from times2saveidx!
    tf(fi,:,:) = abs(as(times2saveidx,[1:50 end-49:end])).^2;
end

% for convenience, compute the difference in power between the two channels
diffmap = squeeze(mean(tf(:,:,1:50),3 )) - squeeze(mean(tf(:,:,51:end),3 ));

%% statistics via permutation testing

% p-value
pval = 0.05;

% convert p-value to Z value
zval = abs(norminv(pval));

% number of permutations
n_permutes = 1000;

% initialize null hypothesis maps
permmaps = zeros(n_permutes,num_frex,length(times2save));

% generate maps under the null hypothesis
for permi = 1:n_permutes
    
    % randomize trials, which also randomly assigns trials to channels
    randorder = randperm(100);
    temp_tf = tf(:,:,randorder);
    
    % compute the "difference" map
    % what is the difference under the null hypothesis?
    permmaps(permi,:,:) = squeeze( mean(temp_tf(:,:,1:50),3) - mean(temp_tf(:,:,51:end),3) );
end

%% show non-corrected thresholded maps

% compute mean and standard deviation maps
mean_h0 = squeeze(mean(permmaps));
std_h0  = squeeze(std(permmaps));

% now threshold real data...
% first Z-score
zmap = (diffmap-mean_h0) ./ std_h0;

% threshold image at p-value, by setting subthreshold values to 0
zmap(abs(zmap)<zval) = 0;


% some plotting...
clim = [0 10000];

figure(1), clf

subplot(221)
imagesc(times2save,frex,diffmap);
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','nor')
title('TF map of real power values')

subplot(222)
imagesc(times2save,frex,zmap);
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-10 10],'xlim',xlim,'ydir','no')
title('Thresholded TF map of Z-values')

subplot(223)
imagesc(times2save,frex,diffmap);
hold on
contour(times2save,frex,logical(zmap),1,'linecolor','k');
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')
title('Power values and outlined significance regions')
colormap jet 

%% corrections for multiple comparisons

% initialize matrices for cluster-based correction
max_cluster_sizes = zeros(1,n_permutes);

% loop through permutations
for permi = 1:n_permutes
    
    % take each permutation map, and transform to Z
    threshimg = squeeze(permmaps(permi,:,:));
    threshimg = (threshimg-mean_h0)./std_h0;
    
    % threshold image at p-value
    threshimg(abs(threshimg)<zval) = 0;
    
    
    % find clusters (need image processing toolbox for this!)
    islands = bwconncomp(threshimg);
    if numel(islands.PixelIdxList)>0
        
        % count sizes of clusters
        tempclustsizes = cellfun(@length,islands.PixelIdxList);
        
        % store size of biggest cluster
        max_cluster_sizes(permi) = max(tempclustsizes);
    end
end

% find cluster threshold (need image processing toolbox for this)
% based on p-value and null hypothesis distribution
cluster_thresh = prctile(max_cluster_sizes,100-(100*pval));

%% plots with multiple comparisons corrections

% now find clusters in the real thresholded zmap
% if they are "too small" set them to zero
islands = bwconncomp(zmap);
for i=1:islands.NumObjects
    % if real clusters are too small, remove them by setting to zero!
    if numel(islands.PixelIdxList{i})<cluster_thresh
        zmap(islands.PixelIdxList{i})=0;
    end
end

% plot tresholded results
figure(2), clf
subplot(221)
imagesc(times2save,frex,diffmap)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('TF power, no thresholding') 
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')

subplot(222)
imagesc(times2save,frex,diffmap)
hold on
contour(times2save,frex,logical(zmap),1,'linecolor','k')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('TF power with contour')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')

subplot(223)
imagesc(times2save,frex,zmap)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('z-map, thresholded')
set(gca,'clim',[-5 5],'xlim',xlim,'ydir','normal')
colormap jet 

%% done.

