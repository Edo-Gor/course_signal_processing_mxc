%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Connectivity
%      VIDEO: All-to-all synchronization and "hubness" (graph theory)
% Instructor: sincxpress.com
%
%%

clear 

% load data
load sampleEEGdata.mat

% apply scalp Laplacian
EEG.data = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);

% pick a frequency to focus on
freq = 6.5; % Hz

% pick a time window for the synchronization analysis
tidx = EEG.times>0 & EEG.times<1000;

%% extract phases from all channels

filtdat = filterFGx( reshape(EEG.data,EEG.nbchan,[]),EEG.srate,freq,3 );
angldat = angle(hilbert( filtdat' ).');
angldat = reshape(angldat,size(EEG.data));

%% all2all synchronization

% initialize
synchmat = zeros(EEG.nbchan);

for chani=1:EEG.nbchan
    for chanj=1:EEG.nbchan
        
        tmpi = squeeze( angldat(chani,tidx,:) );
        tmpj = squeeze( angldat(chanj,tidx,:) );
        
        
        synchmat(chani,chanj) = mean(abs(mean( exp(1i*(tmpi-tmpj)) ,1)),2);
    end
end

%% visualizations

figure(1), clf

clim = [.2 .6];

% the whole matrix
subplot(121)
imagesc(synchmat)
axis square
set(gca,'clim',clim)
xlabel('Channels'), ylabel('Channels')

% topoplot of one channel
chan2plot = 'cz';
subplot(122)
topoplotIndie(synchmat(:,strcmpi({EEG.chanlocs.labels},chan2plot)),EEG.chanlocs,'numcontour',0);
set(gca,'clim',clim)
colormap jet 

%% now for "hubness"

% pick a threshold
thresh = median( nonzeros(triu(synchmat)-eye(EEG.nbchan)) );

% binarize the synchronization image
threshmat = (synchmat-eye(EEG.nbchan)) > thresh;

% compute "hubness"
hubness = sum(threshmat) / (EEG.nbchan-1);


% some visualizations
figure(2), clf
subplot(211)
plot(hubness,'s-','linew',2,'markersize',13,'markerfacecolor','w')
xlabel('Channel number')
ylabel('Hubness (proportion)')

subplot(212)
topoplotIndie(hubness,EEG.chanlocs,'numcontour',0);
colormap jet 

%% exploring "hubness" space

% variable thresholds
thresh = linspace(0,1,37);

hubness = zeros(length(thresh),EEG.nbchan);

% loop over thresholds
for thi=1:length(thresh)
    threshmat = (synchmat-eye(EEG.nbchan)) > thresh(thi);
    hubness(thi,:) = sum(threshmat) / (EEG.nbchan-1);
end

figure(3), clf
h = plot(thresh,hubness);
set(h,'color',ones(3,1)*.4,'marker','o','markerfacecolor','g')
xlabel('"Significance" threshold')
title('Hubness by threshold for all channels')

%% done.
