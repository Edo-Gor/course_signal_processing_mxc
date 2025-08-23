%%
%   COURSE: Neural signal processing and analysis: Zero to hero
%  SECTION: Permutation-based statistics
%  TEACHER: Mike X Cohen, sincxpress.com
%


%% 
% 
%   VIDEO: Permutation testing and shuffling
% 
% 

clear

% number of 'trials' in each condition
n1 = 50;
n2 = 70; % note the trial inbalance!

% create data
data1 = randn(n1,1);
data2 = randn(n2,1) + .3; % note the mean offset!


% step 1: pool the data
alldata = cat(1, data1, data2);

% corresponding labels
truelabels = cat(1, ones(n1,1), 2*ones(n2,1) );

% compute the observed condition difference
true_conddif = mean(alldata(truelabels==2)) - mean(alldata(truelabels==1));

%%% QUESTION: What is the expected mean difference?
%             How does that compare to the empirical mean difference?
%             > it varies a lot

%%


% step 2: shuffle the data once
shuflabels = truelabels( randperm(n1+n2) );


% step 3: compute the mean difference of the shuffled labels
shufconddif = mean(alldata(shuflabels==2)) - mean(alldata(shuflabels==1));

% NOTE that the data have not changed; only the mapping of data to condition label.

%%% QUESTION: How does the shuffled mean difference compare to the observed
%             mean difference?
%             > different
% 
%%% QUESTION: If you re-run this code, how similar are the results?
%             > also quite variable
% 
% 


%% creating a null-hypothesis (H0) distribution

% number of iterations for permutation testing
nIterations = 10000;

permvals = zeros(nIterations,1);


for permi=1:nIterations
    
    % steps 2 and 3 from above
    shuflabels = truelabels(randperm(n1+n2));
    permvals(permi) = mean(alldata(shuflabels==2)) - mean(alldata(shuflabels==1));
end

% show the distribution
figure(1), clf, hold on
histogram(permvals,40)
plot([1 1]*true_conddif,get(gca,'ylim'),'r--','linew',3)
legend({'Shuffled';'True'})
set(gca,'xlim',[-1 1])
xlabel('Mean value')
ylabel('Count')
title(['Permutation test with ' num2str(nIterations) ' iterations'])

%%% QUESTION: How reliable is the distribution?
%             Re-run this code multiple times to observe.
%             Quite reliable
% 
% 
%%% QUESTION: How is the shape of the distribution affected by the number of
%              iterations? Not much overall
%             Try re-running the code with very small and very large number
%              of iterations. What do you think (qualitatively, based on
%              this example) is a reasonable number of iterations?
%              With 100 iteration the shape of the H0 already approximates
%              quite decently a Gaussian 

%% two methods for computing a p-value

% method p_z: compute the normalized distance of the observed p-value from
%             the distribution of H0 p-values

permmean = mean(permvals);
permstd  = std(permvals);

% formula for z-score
zdist = (true_conddif-permmean) / permstd;

% can convert to p-value (1- because it gives back the prob that a value is below -not above- a distribution)
pval = 1 - normcdf(zdist);

subtitle([ 'z_{dist}: ' num2str(zdist) ', p_z value: ' num2str(pval) ])


% method p_c
p_c = sum( permvals>true_conddif ) / nIterations;
subtitle({[ 'z_{dist}: ' num2str(zdist) ', p_z value: ' num2str(pval) ]...
       [ 'p_c value: ' num2str(1-p_c) ]})



%% 
% 
%   VIDEO: Permutation testing in real data
% 
% 


%  The goal here is to determine whether there are statistically
%  significant differences in power between channels 6 and 7.

% load in the v1 data and setup for convolution.

% load in data
clear
load v1_laminar

% useful variables for later...
npnts = size(csd,2);
ntrials = size(csd,3);

chans2use = [ 6 7 ];


% note: too many time points, use post-analysis temporal downsampling!
times2save = -.1:.01:1.14;
times2saveidx = dsearchn(timevec',times2save');


% specify frequency range
min_freq = 10;
max_freq = 90;
num_frex = 35;

% define frequency and range
frex = linspace(min_freq,max_freq,num_frex);


% parameters for complex Morlet wavelets
wavtime  = -2:1/srate:2-1/srate;
half_wav = (length(wavtime)-1)/2;
cycRange = [ 4 10 ];
nCycles  = linspace(cycRange(1),cycRange(2),num_frex);

% FFT parameters
nWave = length(wavtime);
nData = npnts*ntrials;
nConv = nWave+nData-1;

% and create wavelets
cmwX = zeros(num_frex,nConv);
for fi=1:num_frex
    s       = nCycles(fi)/(2*pi*frex(fi)); % frequency-normalized width of Gaussian
    cmw      = exp(1i*2*pi*frex(fi).*wavtime) .* exp( (-wavtime.^2) ./ (2*s^2) );
    tempX     = fft(cmw,nConv);
    cmwX(fi,:) = tempX / max(tempX);
end

%% run convolution to extract tf power. NOTE: You need to keep all the single-trial power

% spectra of data
dataX1 = fft( reshape(csd(chans2use(1),:,:),1,[]) ,nConv);
dataX2 = fft( reshape(csd(chans2use(2),:,:),1,[]) ,nConv);


% initialize output time-frequency data
% notice that here we save all trials
tf = zeros(2,num_frex,length(times2save),ntrials);

for fi=1:num_frex
    
    % run convolution
    as1 = ifft(cmwX(fi,:).*dataX1);
    as1 = as1(half_wav+1:end-half_wav);
    as1 = reshape(as1,npnts,ntrials);
    
    % power on all trials from channel "1"
    % Don't forget to downsample by selecting only the poits of interest
    tf(1,fi,:,:) = abs(as1(times2saveidx,:)).^2; % (note: no averaging b/c keeping all single-trial data!)
    
    
    % run convolution
    as2 = ifft(cmwX(fi,:).*dataX2);
    as2 = as2(half_wav+1:end-half_wav);
    as2 = reshape(as2,npnts,ntrials);
    
    % power on all trials from channel "2"
    tf(2,fi,:,:) = abs(as2(times2saveidx,:)).^2;
end


% for convenience, compute the difference in power between the two channels
diffmap = squeeze(mean(tf(2,:,:,:),4 )) - squeeze(mean(tf(1,:,:,:),4 ));


%% some visualization of the raw power data

clim = [0 20000];
xlim = [-.1 1]; % for plotting

figure(2), clf
subplot(221)
imagesc(times2save,frex,squeeze(mean( tf(1,:,:,:),4 )))
set(gca,'clim',clim,'ydir','n','xlim',xlim)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Channel ' num2str(chans2use(1)) ])
colorbar

subplot(222)
imagesc(times2save,frex,squeeze(mean( tf(2,:,:,:),4 )))
set(gca,'clim',clim,'ydir','n','xlim',xlim)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Channel ' num2str(chans2use(2)) ])
colorbar

subplot(223)
imagesc(times2save,frex,diffmap)
set(gca,'clim',[-mean(clim) mean(clim)],'ydir','n','xlim',xlim)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Difference: channels ' num2str(chans2use(2)) ' - ' num2str(chans2use(1)) ])
colorbar
colormap jet

%% statistics via permutation testing

% p-value
pval = 0.05;

% convert p-value to Z value
% if you don't have the stats toolbox, set zval=1.6449;
zval = abs(norminv(pval));

% number of permutations
n_permutes = 1000;

% initialize null hypothesis maps (h0 tested by shuffling the channels per trials)
permmaps = zeros(n_permutes,num_frex,length(times2save));


% for convenience, tf power maps are concatenated
%   in this matrix, trials 1:ntrials are from channel "1" 
%   and trials ntrials+1:end are from channel "2"
tf3d = cat(3,squeeze(tf(1,:,:,:)),squeeze(tf(2,:,:,:)));



% generate maps under the null hypothesis
for permi = 1:n_permutes
    
    % randomize trials, which also randomly assigns trials to channels
    randorder = randperm(size(tf3d,3));
    temp_tf3d = tf3d(:,:,randorder);
    
    % compute the "difference" map
    % what is the difference under the null hypothesis?
    permmaps(permi,:,:) = squeeze( mean(temp_tf3d(:,:,1:ntrials),3) - mean(temp_tf3d(:,:,ntrials+1:end),3) );
end

%% compute z- and p-values based on normalized distance to H0 distributions (per pixel)

% compute mean and standard deviation maps
mean_h0 = squeeze(mean(permmaps,1));
std_h0  = squeeze(std(permmaps));


%%% now threshold real data...

% first Z-score the real map against the H0 distribution
zmap = (diffmap - mean_h0) ./ std_h0;

% threshold image at p-value, by setting subthreshold values to 0 (notice that we are testing 0.5 on both tails here)
zmap(abs(zmap)<zval) = 0;


%%% now some plotting...

figure(2), clf

subplot(221)
imagesc(times2save,frex,diffmap);
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','nor')
title('TF map of real power values')
colorbar

subplot(222)
imagesc(times2save,frex,zmap);
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-10 10],'xlim',xlim,'ydir','no')
title('Thresholded TF map of Z-values')
colorbar

subplot(223)
imagesc(times2save,frex,diffmap);
hold on
contour(times2save,frex,logical(zmap),1,'linecolor','k');
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')
title('Power values and outlined significance regions')
colorbar
colormap jet


%% 
% 
%    VIDEO: Cluster correction
%
%
%

%%% Note: This code requires having run the code from the previous video.


%  We need to create a distribution of cluster sizes
%  expected under the null hypothesis.


% initialize matrices for cluster-based correction
max_cluster_sizes = zeros(1,n_permutes);


% loop through permutations
for permi = 1:n_permutes
    
    % take each permutation map, and transform to Z
    threshimg = squeeze(permmaps(permi,:,:));
    threshimg = (threshimg - mean_h0) ./ std_h0;
    
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

%% show histograph of maximum cluster sizes

figure(3), clf
hist(max_cluster_sizes,20);
xlabel('Maximum cluster sizes'), ylabel('Number of observations')
title('Expected cluster sizes under the null hypothesis')


% find cluster threshold (need stats toolbox for this line, otherwise find an alternative solution below!)
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
figure(4), clf
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
set(gca,'clim',[-13 13],'xlim',xlim,'ydir','normal')
colormap jet 


%% 
% 
%   VIDEO: Extreme pixel-based multiple comparisons correction
% 
% 
% 

%%% Note: This code requires having run the code from the previous videos.
% Loop again through the permutations, but this time collect the most
%  extreme negative and most extreme positive pixel in the entire map.



% initialize matrices for cluster-based correction
max_vals = zeros(n_permutes,2);


% loop through permutations
for permi = 1:n_permutes
        
    % get extreme values (smallest and largest)
    temp = sort( reshape(permmaps(permi,:,:),1,[] ));
    max_vals(permi,:) = temp([1 end]); % store the largest negative and largest positive values
    
end

% plot the distribution of extreme values
figure(5), clf
histogram(max_vals(:),90)
xlabel('Most extreme H_0 mean difference values')
ylabel('Count')

%% now for thresholding and visualization

pthresh = .025;

% find the lower threshold bound based on percentage
tmp = sort(max_vals(:,1));
thresh_lo = tmp( round(pthresh*n_permutes) );

% repeat for upper threshold bound
tmp = sort(max_vals(:,2));
thresh_hi = tmp( round((1-pthresh)*n_permutes) );


% and plot on top of the histogram
figure(5), hold on
plot([1 1]*thresh_lo,get(gca,'ylim'),'r--','linew',2)
plot([1 1]*thresh_hi,get(gca,'ylim'),'r--','linew',2)


% threshold real data
zmap = diffmap;
% set all pixels between the two thresholds to 0
zmap(zmap>thresh_lo & zmap<thresh_hi) = 0;


figure(6), clf
subplot(221)
imagesc(times2save,frex,diffmap)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('TF power map, no thresholding') 
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','n')

subplot(222), hold on
imagesc(times2save,frex,diffmap)
contour(times2save,frex,logical(zmap),1,'linecolor','k','linewidth',3)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('TF power map with contour')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','normal')

subplot(223)
imagesc(times2save,frex,zmap)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('TF power map, thresholded')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','no')

%%%% optional: different colormap
colormap( bluewhitered(64) )



%% 
% 
%   VIDEO: Extracting features for group analyses
% 
% 
% 


clear

% import data
load groupTFdata.mat
% def n suject for convenience (suubjs - conditions - frex - time)
nSubs = size(tf,1);

% stats parameters
nPerms = 1000;
pVal   = .05; % two-tailed!
sigThresh = norminv(1-pVal/2);


% plot
figure(7), clf
for i=1:2
    subplot(1,2,i)
    contourf(times2save,frex,squeeze(mean(tf(:,i,:,:))),40,'linecolor','none') % average over subjects (1st dim)
    set(gca,'clim',[-1 1]*3)
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title([ 'Condition ' num2str(i) ])
    axis square
    colormap jet
    colorbar
end



% one condition for each subject
figure(8), clf
for i=1:25
    subplot(5,5,i)
    contourf(times2save,frex,squeeze(tf(i,2,:,:)),40,'linecolor','none')
    set(gca,'clim',[-1 1]*3,'xtick',[],'ytick',[])
    axis square
    title([ 's' num2str(i) ])
    colormap jet
    colorbar
end

%% pick a time-frequency region based on condition average

% compute grand average over all subjects and both conditions
grandAve = squeeze(mean(mean(tf,2),1));

% show it
figure(9), clf
contourf(times2save,frex,grandAve,40,'linecolor','none')
set(gca,'clim',[-1 1]*2,'YScale','lo','ytick',ceil(frex(1:8:end)))
axis square
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Grand average over subjects and conditions')
colormap jet
colorbar



%%%% Now, based on visual inspection, draw boxes around two TF windows.
time1 = [285 630]; 
time2 = [180 555];

freq1 = [3 7];
freq2 = [15 30];


% then draw boxes on top
hold on
rectangle('Position',[time1(1) freq1(1) diff(time1) diff(freq1)],'linew',3)
rectangle('Position',[time2(1) freq2(1) diff(time2) diff(freq2)],'linew',3)
text(time1(1),freq1(2),'ROI 1','VerticalAlignment','bottom','FontSize',20)
text(time2(1),freq2(2),'ROI 2','VerticalAlignment','bottom','FontSize',20)


%% Extract data per subject/condition from those ROIs

% convert to indices
time1idx = dsearchn(times2save',time1');
time2idx = dsearchn(times2save',time2');

freq1idx = dsearchn(frex',freq1');
freq2idx = dsearchn(frex',freq2');



%%% extract data
labels = { 'theta cond1';'theta cond2';'beta cond1';'beta cond2' };
data = zeros(4,nSubs);

data(1,:) = squeeze(mean(mean(tf(:,1,freq1idx(1):freq1idx(2),time1idx(1):time1idx(2)),3),4));
data(2,:) = squeeze(mean(mean(tf(:,2,freq1idx(1):freq1idx(2),time1idx(1):time1idx(2)),3),4));
data(3,:) = squeeze(mean(mean(tf(:,1,freq2idx(1):freq2idx(2),freq2idx(1):freq2idx(2)),3),4));
data(4,:) = squeeze(mean(mean(tf(:,2,freq2idx(1):freq2idx(2),freq2idx(1):freq2idx(2)),3),4));



% finally, show a bar plot of the results
figure(10), clf, hold on
bh = bar([1 2 4 5],mean(data,2));
eh = errorbar([1 2 4 5],mean(data,2),std(data,[],2)/sqrt(nSubs-1),'.');
set(gca,'XTickLabel',labels,'xtick',[1 2 4 5])
set(bh,'FaceColor',[.7 .3 .9])
set(eh,'LineWidth',10,'color','k')
ylabel('Power (dB)')

%% done.

