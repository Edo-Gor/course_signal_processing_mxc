%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Connectivity
%      VIDEO: Scalp Laplacian for electrode-level connectivity
% Instructor: sincxpress.com
%
%%

clear 

%% simulate data in two dipoles

% mat file containing EEG, leadfield and channel locations
load emptyEEG


% select dipole location
diploc1 = 109;
diploc2 = 118;

% plot brain dipoles
figure(1), clf, subplot(131)
plot3(lf.GridLoc(:,1), lf.GridLoc(:,2), lf.GridLoc(:,3), 'bo','markerfacecolor','y')
hold on
plot3(lf.GridLoc(diploc1,1), lf.GridLoc(diploc1,2), lf.GridLoc(diploc1,3), 'ks','markerfacecolor','k','markersize',10)
plot3(lf.GridLoc(diploc2,1), lf.GridLoc(diploc2,2), lf.GridLoc(diploc2,3), 'rs','markerfacecolor','r','markersize',10)
rotate3d on, axis square
title('Brain dipole locations')


%% show topomaps before and after Laplacian


% original projections
subplot(232)
topoplotIndie(-lf.Gain(:,1,diploc1), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[-1 1]*40)
title('Dipole 1 projection')

subplot(233)
topoplotIndie(-lf.Gain(:,1,diploc2), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[-1 1]*40)
title('Dipole 2 projection')


% compute laplacians
lapmap1 = laplacian_perrinX(-lf.Gain(:,1,diploc1),[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);
lapmap2 = laplacian_perrinX(-lf.Gain(:,1,diploc2),[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);


% laplacian of projections
subplot(235)
topoplotIndie(lapmap1, EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[-1 1]*400)
title('Laplacian of dipole 1')

subplot(236)
topoplotIndie(lapmap2, EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[-1 1]*400)
title('Laplacian of dipole 2')

colormap jet 

%% compare voltage vs. Laplacian electrode time courses

% don't need so much time...
EEG.times = EEG.times(1:dsearchn(EEG.times',1));
EEG.pnts  = numel(EEG.times);


% generate time series in two dipoles
dipdat = zeros(size(lf.Gain,3),EEG.pnts);
dipdat(diploc1,:) = sin(2*pi*10*EEG.times);
dipdat(diploc2,:) = sin(2*pi*15*EEG.times);
EEG.data = squeeze(lf.Gain(:,1,:))*dipdat;


% also compute laplacian
LAP = EEG;
LAP.data = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);



% Laplacian allows good spatial separation!
plot_simEEG(EEG,27,2);
colormap jet 
plot_simEEG(LAP,27,3);
colormap jet 

% but not everywhere...
plot_simEEG(EEG,31,4);
colormap jet 
plot_simEEG(LAP,31,5);
colormap jet 

%% done.
