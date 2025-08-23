%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Connectivity analyses
%      VIDEO: Power time series correlations
% Instructor: sincxpress.com
%
%%

clear

load v1_laminar.mat

% filter frequencies parameters
minfreq = 10; % Hz
maxfreq = 80; % Hz
numfrex = 30;
frex = linspace(minfreq,maxfreq,numfrex);
widt = linspace(5,10,numfrex);


% time window for analysis
tidx = timevec>0 & timevec<1;

% initialize correlation matrix
cormat = zeros(numfrex,size(csd,1),size(csd,1));

%% do the analysis

% loop over frequencies
for fi=1:numfrex
    
    % amplitude time series
    
    % reshape the data to super-trials
    supdat = reshape(csd,size(csd,1),[]);
    % apply the data
    fdat = filterFGx( supdat,srate,frex(fi),widt(fi) );
    % extract amplitude time series
    powdat = abs(hilbert(fdat')).';
    % reshape back to original size
    powdat = reshape(powdat,size(csd));
    
    % correlation
    
    % loop over trials
    for triali=1:size(powdat,3)
        
        % compute correlations and add to the mix
        cormat(fi,:,:) = squeeze(cormat(fi,:,:)) + ...
                         corr(squeeze(powdat(:,tidx,triali))','type','s');
    end
    
end % end frequencies loop

% scale correlation matrix by trials
cormat = cormat / triali;

%% some plotting

% pick three specific frequencies to plot
frex2plot = [ 12 41 55 ]; % in Hz


figure(1), clf
subplot(211)
plot(frex,mean( reshape(cormat,numfrex,[]).^2 ,2),'ks-','linewid',2,'markerfacecolor','w')
xlabel('Frequency (Hz)')
ylabel('Correlation (R^2)')
title('Sum over all correlations')
ylim([0.075 0.25])

% show a few matrices
for fi=1:3
    
    fidx = dsearchn(frex',frex2plot(fi));
    
    subplot(2,3,3+fi)
    imagesc(squeeze(cormat(fidx,:,:)))
    xlabel('Channels'), ylabel('Channels')
    set(gca,'clim',[-1 1])
    axis square
    title([ ' Correlation matrix at ' num2str(round(frex(fidx))) ' Hz' ])
end

%% prepare an animation

figure(2), clf

% setup the correlation matrix image
subplot(211)
imgh = imagesc(randn(16));
tith = title('asdf');
set(gca,'clim',[-1 1]), axis square
xlabel('Channels'), ylabel('Channels')

% draw the spectral profile
subplot(212), hold on
plot(frex,mean( reshape(cormat,numfrex,[]).^2 ,2),'ks-','linewid',2,'markerfacecolor','w')
lineh = plot([1 1]*frex(1),get(gca,'ylim'),'r--','linew',3);
xlabel('Frequency (Hz)'), ylabel('Correlation (R^2)')
title('Sum over all correlations')
colormap jet


%% animate the animation

for fi=1:numfrex
    
    % update graphics
    set(imgh,'cdata',squeeze(cormat(fi,:,:)))
    set(tith,'string',[ num2str(round(frex(fi),2)) ' Hz' ])
    set(lineh,'xdata',[1 1]*frex(fi))
    
    % update plot
    pause(.5)
end

%% done.
