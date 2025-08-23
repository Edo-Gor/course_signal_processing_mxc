%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Connectivity analyses
%      VIDEO: Project 5-2: Seeded phase vs. power coupling
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
cormat = zeros(2,numfrex,size(csd,1),size(csd,1));

%% do the analysis

% loop over frequencies
for fi=1:numfrex
    
    %%% amplitude time series
    
    % reshape the data to super-trials
    supdat = reshape(csd,size(csd,1),[]);
    % apply the data
    hdat = hilbert( filterFGx( supdat,srate,frex(fi),widt(fi) )' ).';
    hdat = reshape(hdat,size(csd));
    hdat = hdat(:,tidx,:);
    
    % extract amplitude and phase angle time series
    powdat = abs(hdat);
    phsdat = angle(hdat);
    
    %%% correlation
    
    % loop over trials
    for triali=1:size(powdat,3)
        
        % compute amplitude time series correlations
        cormat(1,fi,:,:) = squeeze(cormat(1,fi,:,:)) + corr(squeeze(powdat(:,:,triali))','type','s');
        
        % compute phase synchronization
        for chani=1:size(csd,1)
            cormat(2,fi,chani,:) = squeeze(cormat(2,fi,chani,:)) + abs(mean(exp(1i*bsxfun(@minus,phsdat(:,:,triali),phsdat(chani,:,triali))),2));
        end
    end
    
end % end frequencies loop

% scale correlation matrix by trials
cormat = cormat / triali;

%% some plotting

% pick three specific frequencies to plot
frex2plot = [ 12 41 55 ]; % in Hz


figure(1), clf
subplot(311), hold on
plot(frex,mean( reshape(cormat(1,:,:,:),numfrex,[]).^2 ,2),'ks-','linewid',2,'markerfacecolor','w')
plot(frex,mean( reshape(cormat(2,:,:,:),numfrex,[]).^2 ,2),'ro-','linewid',2,'markerfacecolor','w')
xlabel('Frequency (Hz)')
ylabel('Correlation or synchronization')
title('Sum over all correlations')
legend({'Power time series';'Phase synchronization'})

% show a few matrices
for fi=1:3
    
    fidx = dsearchn(frex',frex2plot(fi));
    
    % power correlations
    subplot(3,3,3+fi)
    imagesc(squeeze(cormat(1,fidx,:,:)))
    xlabel('Channels'), ylabel('Channels')
    set(gca,'clim',[-1 1])
    axis square
    title([ 'Power corrs.: ' num2str(round(frex(fidx))) ' Hz' ])
    
    
    % phase synchronization
    subplot(3,3,6+fi)
    imagesc(squeeze(cormat(2,fidx,:,:)))
    xlabel('Channels'), ylabel('Channels')
    set(gca,'clim',[0 1])
    axis square
    title([ 'Phase synch.: ' num2str(round(frex(fidx))) ' Hz' ])
    colormap jet
end

%% show "topography"

% select plotting parameters
freq2plot = 50;
seedchan  = [2 6 15];

% now for the plotting
figure(2), clf
subplot(131)
imagesc(squeeze(cormat(:,dsearchn(frex',freq2plot),seedchan(1),:))')
title({'Seeded connectivity:';['channel ' num2str(seedchan(1)) ' at ' num2str(freq2plot) ' Hz' ]})
set(gca,'xtick',1:2,'xticklabel',{'Power';'Phase'})
ylabel('Channel')
axis image
set(gca,'clim',[-1 1])
colorbar

subplot(132)
imagesc(squeeze(cormat(:,dsearchn(frex',freq2plot),seedchan(2),:))')
title({'Seeded connectivity:';['channel ' num2str(seedchan(2)) ' at ' num2str(freq2plot) ' Hz' ]})
set(gca,'xtick',1:2,'xticklabel',{'Power';'Phase'})
ylabel('Channel')
axis image
set(gca,'clim',[-1 1])
colorbar

subplot(133)
imagesc(squeeze(cormat(:,dsearchn(frex',freq2plot),seedchan(3),:))')
title({'Seeded connectivity:';['channel ' num2str(seedchan(3)) ' at ' num2str(freq2plot) ' Hz' ]})
set(gca,'xtick',1:2,'xticklabel',{'Power';'Phase'})
ylabel('Channel')
axis image
set(gca,'clim',[-1 1])
colorbar

colormap jet

%% done.
