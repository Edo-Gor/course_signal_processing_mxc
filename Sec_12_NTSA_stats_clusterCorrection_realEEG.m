%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Inferrential statistics
%      VIDEO: Permutation testing for cluster correction in real EEG data
% Instructor: sincxpress.com
%
%%

clear

% import data
load groupTFdata.mat
% convenience
nSubs = size(tf,1);

% stats parameters
nPerms = 1000;
pVal   = .05; % two-tailed!
sigThresh = norminv(1-pVal/2);


% plot
figure(1), clf
for i=1:2
    subplot(1,2,i)
    contourf(times2save,frex,squeeze(mean(tf(:,i,:,:))),40,'linecolor','none')
    set(gca,'clim',[-1 1]*3)
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title([ 'Condition ' num2str(i) ])
    axis square
    colormap jet
end



% one condition for each subject
figure(2), clf
for i=1:25
    subplot(5,5,i)
    contourf(times2save,frex,squeeze(tf(i,2,:,:)),40,'linecolor','none')
    set(gca,'clim',[-1 1]*3,'xtick',[],'ytick',[])
    axis square
    title([ 's' num2str(i) ])
    colormap jet
end

%% permutation testing

% compute true differences
tfdiff = squeeze(diff(tf,[],2));
tfFake = zeros(size(tfdiff));


% initialize
permutedDiffs = zeros(nPerms,size(tf,3),size(tf,4));

for permi = 1:nPerms
    
    % random sequence of +/-1
    fakeconds = sign(randn(nSubs,1));
    
    % randomize sign of difference map
    for si=1:nSubs
        tfFake(si,:,:) = fakeconds(si)*tfdiff(si,:,:);
    end
    
    % compute and store difference time series
    permutedDiffs(permi,:,:) = mean(tfFake);
end

%% pixel-based statistical thresholding

[zmap,zthreshMap] = deal( squeeze( (mean(tfdiff)-mean(permutedDiffs,1)) ./ std(permutedDiffs,[],1) ) );
zthreshMap(abs(zmap)<sigThresh) = 0;


clim = [-1 1]*5;

figure(3), clf
subplot(221)
contourf(times2save,frex,squeeze(mean(tfdiff)),40,'linecolor','none')
set(gca,'clim',[-1 1]*2)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('"Raw" differences')
axis square

subplot(222)
contourf(times2save,frex,zmap,40,'linecolor','none')
set(gca,'clim',clim)
title('Z-map')
axis square

subplot(223)
contourf(times2save,frex,zthreshMap,40,'linecolor','none')
set(gca,'clim',clim)
title('Thresholded z-map')
axis square
colormap jet 

%% find cluster sizes under the null hypothesis

% initialize cluster sizes from permutation
clustsizes = zeros(nPerms,1);

for permi=1:nPerms
    
    % compute z-score difference
    zdiffFake = squeeze( (permutedDiffs(permi,:,:)-mean(permutedDiffs)) ./ std(permutedDiffs) );
    
    % threshold
    zdiffFake( abs(zdiffFake)<sigThresh ) = 0;
    
    % identify clusters
    islands = bwconncomp( logical(zdiffFake) );
    
    % find cluster sizes
    clustNs = cellfun(@length,islands.PixelIdxList);
    clustsizes(permi) = max(clustNs);
end

% compute cluster threshold
clustthresh = prctile(clustsizes,100-pVal*100);


% show distribution of cluster sizes
figure(4), clf, hold on
histogram(clustsizes)
plot([1 1]*clustthresh,get(gca,'ylim'),'r--','linew',3)
xlabel('Cluster size (time points)'), ylabel('Count')


%% remove small clusters from real thresholded data

% find islands
islands = bwconncomp( logical(zthreshMap) );

% find and remove any subthreshold islands
for ii=1:islands.NumObjects
    if numel(islands.PixelIdxList{ii})<clustthresh
        zthreshMap(islands.PixelIdxList{ii}) = 0;
    end
end

% now plot that
figure(3)
subplot(224), hold on
contourf(times2save,frex,zmap,40,'linecolor','none')
contour (times2save,frex,logical(zthreshMap),1,'linecolor','k')
title('Statistical results, cluster-corrected')
set(gca,'clim',clim*2)
axis square
colormap jet 


%% done!
