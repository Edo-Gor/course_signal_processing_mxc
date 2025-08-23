%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Simulating EEG data
%      VIDEO: Project 2: dipole-level EEG data
% Instructor: sincxpress.com
%
%%

%% 

% mat file containing EEG, leadfield and channel locations
load emptyEEG

% select dipole location (more-or-less random)
diploc = 109;

% plot brain dipoles
figure(1), clf, subplot(121)
plot3(lf.GridLoc(:,1), lf.GridLoc(:,2), lf.GridLoc(:,3), 'bo','markerfacecolor','y')
hold on
plot3(lf.GridLoc(diploc,1), lf.GridLoc(diploc,2), lf.GridLoc(diploc,3), 'rs','markerfacecolor','k','markersize',10)
rotate3d on, axis square
title('Brain dipole locations')


% Each dipole can be projected onto the scalp using the forward model. 
% The code below shows this projection from one dipole.
subplot(122)
topoplotIndie(-lf.Gain(:,1,diploc), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[-1 1]*40)
title('Signal dipole projection')

%% add signal to one dipole and project to scalp

% reduce data size a bit
EEG.pnts  = 2000;
EEG.times = (0:EEG.pnts-1)/EEG.srate;

% initialize all dipole data
dipole_data = zeros(size(lf.Gain,3),EEG.pnts);

% add signal to one dipole
dipole_data(diploc,:) = sin(2*pi*10*EEG.times);

% now project dipole data to scalp electrodes
EEG.data = squeeze(lf.Gain(:,1,:))*dipole_data;

% plot the data
plot_simEEG(EEG,31,2);

%% now for the projects!

%%%% IMPORTANT! Check internal consistency with existing EEG structure!!!

EEG


%% 1) pure sine wave with amplitude explorations

EEG.trials = 40;

% dipole amplitude magnitude
ampl = 1/10;

% initialize all dipole data
dipole_data = zeros(size(lf.Gain,3),EEG.pnts);
dipole_data(diploc,:) = ampl*sin(2*pi*10*EEG.times);

% compute one dipole
EGG.data = squeeze(lf.Gain(:,1,:))*dipole_data;

% repeat that for N trials
% you can replace the loop with the following line:
% EEG.data = repmat(EEG.data,[1 1 EEG.trials]);
for triali = 1:EEG.trials
    EEG.data(:,:,triali) = squeeze(lf.Gain(:,1,:))*dipole_data;
end

% plot the data
plot_simEEG(EEG,31,3);

%%% Question: What is the smallest amplitude of dipole signal that still
%             elicits a scalp-level response? 
%             Your basically never getting there becouse all the other
%             electrodes are set to 0

%% 2) sine wave with noise

%%% Question: Given amplitude=1 of dipole signal, what standard deviation of noise
%             at all other dipoles overpowers the signal (qualitatively)?
%             From ampln = 1/2 it starts to be difficult but al least the
%             spectrum is meaningful, from 1 on, it's gone

% dipole amplitude magnitude
ampl = 1;
% noise ampliude magnitude
ampln = 1/10;

% initialize all dipole data
dipole_data = zeros(size(lf.Gain,3),EEG.pnts);

% repeat that for N trials
for triali = 1:EEG.trials
    
    % dipole data (inside so to loop over randn)
    dipole_data(:,:) = ampln*randn(size(dipole_data));
    dipole_data(diploc,:) = ampl*sin(2*pi*10*EEG.times);
    
    EEG.data(:,:,triali) = squeeze(lf.Gain(:,1,:))*dipole_data;
end

% plot the data
plot_simEEG(EEG,31,4);


%% 3) Non-oscillatory transient in one dipole, noise in all other dipoles

% noise ampliude magnitude
ampln = 1;

% initialize all dipole data
dipole_data = zeros(size(lf.Gain,3),EEG.pnts);

% set the gaussian
peaktime = 1; % seconds
width = .10;
ampl = 70;
gaus = ampl * exp( -(EEG.times-peaktime+rand/5).^2 / (2*width^2) );

% repeat that for N trials
for triali = 1:EEG.trials
     
    % dipole data (inside so to loop over randn)
    dipole_data(:,:) = ampln*randn(size(dipole_data));
    dipole_data(diploc,:) = gaus;
    
    EEG.data(:,:,triali) = squeeze(lf.Gain(:,1,:))*dipole_data;
end

plot_simEEG(EEG,31,5);

%% 4) Non-stationary oscillation in one dipole, transient oscillation in another dipole, noise in all dipoles

%%% first pick two dipoles
diploc1 = 110;
diploc2 = 510;

% plot brain dipoles
figure(6), clf, subplot(2,2,[1,3])
plot3(lf.GridLoc(:,1), lf.GridLoc(:,2), lf.GridLoc(:,3), 'bo','markerfacecolor','y')
hold on
plot3(lf.GridLoc(diploc1,1), lf.GridLoc(diploc1,2), lf.GridLoc(diploc1,3), 'rs','markerfacecolor','k','markersize',10)
plot3(lf.GridLoc(diploc,1), lf.GridLoc(diploc,2), lf.GridLoc(diploc,3), 'ks','markerfacecolor','r','markersize',10)
rotate3d on, axis square
title('Brain dipole locations')


% Each dipole can be projected onto the scalp using the forward model. 
% The code below shows this projection from one dipole.
subplot(222)
topoplotIndie(-lf.Gain(:,1,diploc1), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[-1 1]*40)
title('Signal dipole1 projection')

subplot(224)
topoplotIndie(-lf.Gain(:,1,diploc2), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
set(gca,'clim',[-1 1]*40)
title('Signal dipole2 projection')



%%% then do the simulation

% amplitude magnitude for aoverall noise
ampln = 1/5;

% initialize all dipole data
dipole_data = zeros(size(lf.Gain,3),EEG.pnts);

% set the non-stationary oscillation for dipole 110


% set the gaussian for dipole 510
peaktime = 1; % seconds
width = .10;
sinefreq = 7;
ampl = 1;
gaus = ampl * exp( -(EEG.times-peaktime+rand/5).^2 / (2*width^2) );

% repeat that for N trials
for triali = 1:EEG.trials
     
    % dipole data (inside so to loop over randn)
    dipole_data(:,:) = ampln*randn(size(dipole_data));
    % dipole 1 non-stationary oscillation (range 5 - 10 Hz); scaled by 2
    % to highlight it
    freqmod = 5 + 5*interp1(rand(1,10),linspace(1,10,EEG.pnts));
    dipole_data(diploc1,:) = 2*sin( 2*pi * ((EEG.times + cumsum(freqmod))/EEG.srate) );
    % dipole 1 transient oscillation (wave by gauss taper) 
    dipole_data(diploc2,:) = sin(2*pi*sinefreq*EEG.times + pi*rand) .* gaus;
    
    EEG.data(:,:,triali) = squeeze(lf.Gain(:,1,:))*dipole_data;
end

% plot the data
plot_simEEG(EEG,56,7);
plot_simEEG(EEG,31,8);

%% 
