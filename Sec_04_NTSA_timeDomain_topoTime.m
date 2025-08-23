%%
%     COURSE: Solved problems in neural time series analysis
%    SECTION: Time-domain analyses
%      VIDEO: Topography time series
% Instructor: sincxpress.com
%
%%

load sampleEEGdata.mat

% time points for topographies
times2plot = -200:50:800; % in ms

% convert to indices (search for the EEG timepoint whose index is closest
% to our plot scale)
tidx = zeros(size(times2plot));
for i=1:length(tidx)
    [~,tidx(i)] = min(abs( EEG.times-times2plot(i) ));
end

% or (notice that the search is actually faster with the loop):
tidx = dsearchn(EEG.times',times2plot');

%% topoplot time series at exact time point

figure(1), clf

% optional: redo with average reference


% define subplot geometry (how many rows and cols are needed)
subgeomR = ceil(sqrt(length(tidx)));
subgeomC = ceil(length(tidx)/subgeomR);

% compute ERP
erp = mean(EEG.data,3);

for i=1:length(tidx)
    subplot( subgeomR,subgeomC,i )
    topoplotIndie(erp(:,tidx(i)),EEG.chanlocs );
    set(gca,'clim',[-1 1]*10)  % otherwise each map has its own colormap
    title([ num2str(times2plot(i)) ' ms' ])
end

%% topoplot time series at average around time point
% in order to get less noisy pictures, more robust

% window size
twin = 10; % in ms; half of window

% convert to indices
twinidx = round(twin/(1000/EEG.srate));

figure(2), clf
for i=1:length(tidx)
    subplot( subgeomR,subgeomC,i )
    
    % time points to average together
    times2ave = tidx(i)-twinidx : tidx(i)+twinidx;
    
    % draw the topomap
    topoplotIndie( mean(mean(EEG.data(:,times2ave,:),3),2),EEG.chanlocs,'electrodes','off','numcontour',0 );
    set(gca,'clim',[-1 1]*10)
    title([ num2str(times2plot(i)) ' ms' ])
end

%% topoplot time series at exact time point
% optional: redo with average reference

figure(3), clf

% average ref
EEG.cardata = bsxfun(@minus,EEG.data,mean(EEG.data,1));

% define subplot geometry (how many rows and cols are needed)
subgeomR = ceil(sqrt(length(tidx)));
subgeomC = ceil(length(tidx)/subgeomR);

% compute ERP
erp = mean(EEG.cardata,3);

for i=1:length(tidx)
    subplot( subgeomR,subgeomC,i )
    topoplotIndie(erp(:,tidx(i)),EEG.chanlocs );
    set(gca,'clim',[-1 1]*10)  % otherwise each map has its own colormap
    title([ num2str(times2plot(i)) ' ms' ])
end


%% done.
