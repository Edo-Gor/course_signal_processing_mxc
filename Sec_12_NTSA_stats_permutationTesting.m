%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Inferrential statistics
%      VIDEO: Permutation testing in simulated time series
% Instructor: sincxpress.com
%
%%

clear

% generic simulation parameters
npnts = 2501;
time  = linspace(-1,3,npnts)'; % force column vector
srate = 1./mean(diff(time));   % Hz (seems weird like this but helpful to filter out small irregularities in the emprirical sampling rate)
ntrials = 100;  % per condition!
nPerms  = 1000; % number of iterations in permutation testing
pval    = .05;  % p-value threshold

% signal-specific parameters
amplit1 = 4;
amplit2 = 3;
noise1  = 2;
noise2  = 2;

% create two "pure" signals
pure1 = amplit1   * exp( -4*log(2)*(time-1).^2 / .25 ) + ...
        amplit1/2 * exp( -4*log(2)*(time-2).^2 / .05 );

pure2 = amplit2   * exp( -4*log(2)*(time-1).^2 / .25 ) + ...
        amplit2/2 * exp( -4*log(2)*(time-2).^2 / .05 );


% create dataset: pure signal plus random noise (take outer product to avoid a loop)
data1 = pure1*ones(1,ntrials) + randn(npnts,ntrials)*noise1;
data2 = pure2*ones(1,ntrials) + randn(npnts,ntrials)*noise2;


% plot
figure(1), clf
subplot(311)
plot(time,pure1, time,pure2)
ylabel('Amplitude (a.u.)')
legend({'signal1';'signal2'})
title('Pure signals')

subplot(312)
plot(time,mean(data1,2), time,mean(data2,2))
ylabel('Amplitude (a.u.)')
legend({'data1';'data2'})
title([ 'Signals + noise (average of ' num2str(ntrials) ' trials)' ])

%% permutation testing

% put data together into one matrix with 2N trials
data3d = cat(2,data1,data2);

% generate true condition labels
condlabels = (1:ntrials*2)>ntrials;

% initialize
permutedDiffs = zeros(npnts,nPerms);

for permi = 1:nPerms
    
    % shuffle condition label vector
    fakeconds = condlabels( randperm(2*ntrials) );
    
    % compute and store difference time series
    mean1 = mean( data3d(:,fakeconds==0),2 );
    mean2 = mean( data3d(:,fakeconds==1),2 );
    permutedDiffs(:,permi) = mean1-mean2;
end


% compute z-score difference
obsdiff = mean(data3d(:,condlabels==0),2) - mean(data3d(:,condlabels==1),2);
zdiff   = (obsdiff-mean(permutedDiffs,2)) ./ std(permutedDiffs,[],2);

% statistically threshold the final result
sigThresh = norminv(1-pval/2); % note: two-tailed!
zthresh = zdiff;
zthresh( abs(zthresh)<sigThresh ) = 0;

%% plot 

subplot(313), hold on
plot(time,zdiff,'k')
plot(time(logical(zthresh)), zthresh(find(zthresh)),'yo','markerfacecolor','r')
xlabel('Time (s)'), ylabel('Z value')
title('Statistical results')
legend({'Difference';[ 'p<' num2str(pval) ', 2-tailed' ]})

%% 
