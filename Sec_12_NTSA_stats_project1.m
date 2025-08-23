%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Inferrential statistics
%      VIDEO: Project 7-1: Effects of noise smoothness on cluster correction
% Instructor: sincxpress.com
%
%%

clear

% generic simulation parameters
npnts   = 2501;
time    = linspace(-1,3,npnts);
srate   = 1./mean(diff(time)); % Hz
ntrials =  100; % per condition!
nPerms  = 1000; % number of iterations in permutation testing

% convert p-value to z-value (note: use 1.96 for p<.05 2-tailed without stats-toolbox)
pval      =  .05; % p-value threshold
sigThresh = sqrt(2) * erfinv(1 - pval);

% signal-specific parameters
amplit1 = 4;
amplit2 = 3;
fwhm1   = .5;
fwhm2   = .3;

% create two "pure" signals
pure1 = amplit1   * exp( -4*log(2)*(time-1).^2 / fwhm1^2 ) + ...
        amplit1/1 * exp( -4*log(2)*(time-2).^2 / fwhm2^2 );

pure2 = amplit2   * exp( -4*log(2)*(time-1).^2 / fwhm1^2 ) + ...
        amplit2/2 * exp( -4*log(2)*(time-2).^2 / fwhm2^2 );


%% create dataset with noise

% noise parameters
noisestd = 2;
peakfreq = 4;
fwhm     = 5;

%%% create frequency-domain Gaussian
hz = linspace(0,srate,npnts);
s  = fwhm*(2*pi-1)/(4*pi); % normalized width
x  = hz-peakfreq;          % shifted frequencies
fg = exp(-.5*(x/s).^2);    % gaussian



% create dataset: pure signal plus noise
[data1,data2] = deal( zeros(npnts,ntrials) );
for triali=1:ntrials
    
    % data1
    fc = rand(1,npnts) .* exp(1i*2*pi*rand(1,npnts));
    noise = real( ifft(fc.*fg) );
    noise = noise*(noisestd/std(noise)); % adjust STD
    data1(:,triali) = pure1 + noise;
    
    % data2
    fc = rand(1,npnts) .* exp(1i*2*pi*rand(1,npnts));
    noise = real( ifft(fc.*fg) );
    noise = noise*(noisestd/std(noise)); % adjust STD
    data2(:,triali) = pure2 + noise;
end


% plot
figure(1), clf
subplot(211)
plot(time,pure1, time,pure2)
ylabel('Amplitude (a.u.)')
legend({'signal1';'signal2'})
title('Pure signals')

subplot(212)
plot(time,mean(data1,2), time,mean(data2,2))
ylabel('Amplitude (a.u.)')
legend({'data1';'data2'})
title([ 'Signals + noise (average of ' num2str(ntrials) ' trials)' ])


%% permutation testing

% put data together into one matrix with 2N trials
data3d = cat(2,data1,data2);

% generate true condition labels
condlabels = (1:ntrials*2)>ntrials;

permutedDiffs = zeros(npnts,nPerms);
for permi = 1:nPerms
    
    % shuffle condition label vector
    fakeconds = condlabels( randperm(2*ntrials) );
    
    % compute and store difference time series
    mean1 = mean( data3d(:,fakeconds==0),2 );
    mean2 = mean( data3d(:,fakeconds==1),2 );
    permutedDiffs(:,permi) = mean1-mean2;
    
end

%% find cluster sizes under the null hypothesis


% initialize cluster sizes from permutation
clustsizes = zeros(nPerms,1);

for permi=1:nPerms
    
    % compute z-score difference
    zdiffFake = (permutedDiffs(:,permi)-mean(permutedDiffs,2)) ./ std(permutedDiffs,[],2);
    
    % threshold
    zdiffFake( abs(zdiffFake)<sigThresh ) = 0;
    
    % identify clusters
    islands = bwconncomp( logical(zdiffFake) );
    
    % find cluster sizes
    clustNs = cellfun(@length,islands.PixelIdxList);
    
    if ~isempty(clustNs)
        clustsizes(permi) = max(clustNs);
    end
end

% compute cluster threshold
clustthresh = prctile(clustsizes,100-pval*100);


% show distribution of cluster sizes
figure(2), clf
histogram(clustsizes)
xlabel('Cluster size (time points)'), ylabel('Count')

%% remove small clusters from real thresholded data

% compute z-score difference
obsdiff = mean(data3d(:,condlabels==0),2) - mean(data3d(:,condlabels==1),2);
zdiff   = (obsdiff-mean(permutedDiffs,2)) ./ std(permutedDiffs,[],2);

% recompute thresholded time series
zthresh = zdiff;
zthresh( abs(zthresh)<sigThresh ) = 0;

% plot that
figure(3), clf
subplot(311), hold on
plot(time,zdiff,'k')
plot(time(logical(zthresh)), zthresh(logical(zthresh)),'yo','markerfacecolor','r')
ylabel('Z value')
title('Statistical results, uncorrected')
legend({'Difference';[ 'p<' num2str(pval) ', 2-tailed' ]})




% find islands
islands = bwconncomp( logical(zthresh) );

% find and remove any subthreshold islands
for ii=1:islands.NumObjects
    if numel(islands.PixelIdxList{ii})<clustthresh
        zthresh(islands.PixelIdxList{ii}) = 0;
    end
end

% now plot that
subplot(312), hold on
plot(time,zdiff,'k')
plot(time(logical(zthresh)), zthresh(logical(zthresh)),'yo','markerfacecolor','r')
ylabel('Z value')
title('Statistical results, cluster-corrected')
legend({'Difference';[ 'p<' num2str(pval) ', 2-tailed' ]})



% now plot that nicer
subplot(313), cla, hold on
plot(time,zdiff,'k')

islands = bwconncomp( logical(zthresh) );
for ii=1:islands.NumObjects
    iidx = islands.PixelIdxList{ii};
    patch([ time(iidx) time(iidx(end-1:-1:1)) ],[zdiff(iidx)' zeros(1,numel(iidx)-1)],'r','edgecolor','none');
end

xlabel('Time (s)'), ylabel('Z value')
title('Statistical results, corrected')
legend({'Difference';[ 'p<' num2str(pval) ', 2-tailed' ]})

%% done!
