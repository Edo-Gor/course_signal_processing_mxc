%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Inferrential statistics
%      VIDEO: Permutation testing for one variable and two groups
% Instructor: sincxpress.com
%
%%

clear

%% generate data

% dataset parameters
mean1 = 500; % ms
mean2 = 480; % ms
stdev =  60; % ms

nTrials = 100; % per condition!

% first column is condition number
data(:,1) = double( (1:2*nTrials)>nTrials );

% second column is RTs (random numbers with std=stdev and mean=mean1/2)
data(1:nTrials,2)     = mean1 + randn(nTrials,1) * stdev;
data(nTrials+1:end,2) = mean2 + randn(nTrials,1) * stdev;

%% visualize data distributions and their difference

figure(1), clf
subplot(5,1,2:4), hold on
[y1,x1] = hist(data(data(:,1)==0,2),20); % histogram of condition=0
[y2,x2] = hist(data(data(:,1)==1,2),20); % histogram of condition=1

% plot the distributions
h = plot(x1,y1, x2,y2, 'linew',3);
set(h,'linewidth',3)
set(h(1),'color','b')
set(h(2),'color','r')

% plot the means
plot([1 1]*mean(data(data(:,1)==0,2)),get(gca,'ylim'),'b--')
plot([1 1]*mean(data(data(:,1)==1,2)),get(gca,'ylim'),'r--')

% some labels and legends
xlabel('RTs (ms)'), ylabel('Count')
legend('C1 dist.','C2 dist.','C1 mean','C2 mean')

%% permutation test

nPermutes = 1000;
permMeans = zeros(nPermutes,1);

trueCondLabs = data(:,1);

for permi=1:nPermutes
    
    % shuffle condition labels
    fakeCondLabs = trueCondLabs(randperm(length(trueCondLabs)));
    
    % compute fake means
    fakemean1 = mean(data(fakeCondLabs==0,2));
    fakemean2 = mean(data(fakeCondLabs==1,2));
    
    % permutation is their difference
    permMeans(permi) = fakemean1 - fakemean2;
end

% compute observed mean difference
realmean1 = mean(data(data(:,1)==0,2));
realmean2 = mean(data(data(:,1)==1,2));
realdiff  = realmean1 - realmean2;

% compare real value against H0 distribution
zval = (realdiff - mean(permMeans)) / std(permMeans);  % hint: z = (obs-center)/std
pval = normcdf(-abs(zval));

% compute parametric t-statistic
[~,p,~,t] = ttest2(data(data(:,1)==0,2),data(data(:,1)==1,2));
% convert t into z
t2p = tpdf(t.tstat,t.df);
t2z = sign(zval) * norminv(1-t2p);


%% show null-hypothesis distribution

figure(2), clf, hold on

% permuted distribution
hist(permMeans,50)
xlabel('Mean difference (ms)')
ylabel('Count')

% observed difference
plot([1 1]*realdiff,get(gca,'ylim'),'r--','linew',4)

title({[ 'Z_{perm}=' num2str(round(zval,2)) , ', p_{perm}=' num2str(round(pval,3)) ],...
       [ 'Z_{para}=' num2str(round(t2z,2))  , ', p_{para}=' num2str(round(t2p,3))  ] })

%% done.
