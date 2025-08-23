%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Connectivity
%      VIDEO: Phase lag index
% Instructor: sincxpress.com
%
%%

clear 

%% basic usage

N = 100;

phases = pi*rand(N,1) + pi/2;

% compute PLI
euler  = exp(1i* phases );
signIm = sign(imag( euler ));
pli    = abs(mean(signIm));

% compute ISPC
ispc   = abs(mean(euler));


%%% plotting
figure(1), clf
h = polarplot([zeros(N,1) phases]',[zeros(N,1) ones(N,1)]');
set(h,'color',ones(1,3)*.4)
hold on
h = polarplot([0 angle(mean(euler))],[0 abs(mean(euler))]);
set(h,'color','k','linew',6)

th = title([ 'PLI: ' num2str(pli) ', ISPC: ' num2str(ispc) ]);
set(th,'FontSize',15)

%% prepare simulation showing varying phase lags

% basic stuff
pnts = 300;
time = (0:pnts-1)/pnts;

% phase angle differences
phases = linspace(0,2*pi,100);

% generate first signal (static)
signal1 = sin( linspace(0,6*pi,300) ) + randn(1,pnts)/5;
angles1 = angle(hilbert(signal1));

% initialize results vectors
[plis,ispcs] = deal( zeros(1,length(phases)) );



%%% setup figure
figure(2), clf
subplot(211)
sineh = plot(time,signal1, time,signal1);
set(gca,'ylim',[-1 1]*1.5)
set(sineh,'linew',2)
xlabel('Time (a.u.)')
legend({'Signal 1';'Signal 2'})

% polar plot for phase angle differences
subplot(223)
polarh = polarplot([zeros(1,pnts); rand(1,pnts)],([0 1]'*ones(1,pnts)));
hold on
meanVh = polarplot([0 1],[0 1]);
set(polarh,'color',ones(3,1)*.7)
set(meanVh,'color','k','linew',3)

% connectivity measures
subplot(224)
connh = plot(phases,plis, phases,ispcs);
set(gca,'xlim',phases([1 end]))
set(connh,'linew',2)
xlabel('Phases (rad.)'), ylabel('Connectivity strength')
legend({'PLI';'ISPC'})


%% now for the animation

% loop over phases
for phi=1:length(phases)
    
    % generate second signal
    signal2 = sin( linspace(0,6*pi,300) + phases(phi) ) + randn(1,pnts)/5;
    angles2 = angle(hilbert(signal2));
    
    % Eulerized phase differences
    phaseDiffs = exp(1i* (angles1-angles2) );
    
    % ISPC
    ispcs(phi) = abs(mean( phaseDiffs ));
    % PLI
    plis(phi)  = abs(mean(sign(imag( phaseDiffs ))));
    
    
    
    %%% adjust figure
    % signal 2
    set(sineh(2),'ydata',signal2)
    
    % polar plot
    for ii=1:pnts, set(polarh(ii),'ThetaData',[0 angle(phaseDiffs(ii))]), end
    set(meanVh,'ThetaData',[0 angle(mean(phaseDiffs))], 'RData',[0 abs(mean(phaseDiffs))]);
    
    % connectivity plot
    set(connh,'xdata',phases(1:phi))
    set(connh(1),'ydata',plis(1:phi))
    set(connh(2),'ydata',ispcs(1:phi))
    
    % update
    pause(.1)
end

%% done.
