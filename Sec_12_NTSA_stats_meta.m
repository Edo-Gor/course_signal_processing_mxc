%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Inferrential statistics
%      VIDEO: Meta-permutation test for increased stability
% Instructor: sincxpress.com
%
%%

%% generate data

% dataset parameters
mean1 = 500; % ms
mean2 = 485; % ms
stdev =  50; % ms

nTrials = 100; % per condition!

% first column is condition number
data(:,1) = double( (1:2*nTrials)>nTrials );

% second column is RTs
data(1:nTrials,2)     = stdev*randn(nTrials,1) + mean1;
data(nTrials+1:end,2) = stdev*randn(nTrials,1) + mean2;

%% permutation test

nPermRuns = 100;
nPermutes = 1000;

permMeans = zeros(nPermutes,nPermRuns);


for metai=1:nPermRuns
    
    for permi=1:nPermutes
        
        % shuffle condition labels
        fakeCondLabs = randn(nTrials*2,1)>0;
        
        % compute fake means
        fakemean1 = mean(data(fakeCondLabs==0,2));
        fakemean2 = mean(data(fakeCondLabs==1,2));
        
        % permutation is their difference
        permMeans(permi,metai) = fakemean1 - fakemean2;
    end
    
end

% compute observed mean difference
realmean1 = mean(data(data(:,1)==0,2));
realmean2 = mean(data(data(:,1)==1,2));
realdiff  = realmean1 - realmean2;

% compare real value against H0 distribution (note that here you deal with vectors)
zval = (realdiff-mean(permMeans)) ./ std(permMeans);
pval = normcdf(-abs(zval));

%% show all null-hypothesis distributions

figure(1), clf, hold on

[x,y] = deal( zeros( 50,nPermRuns ) );
for i=1:nPermRuns
    [y(:,i),x(:,i)] = hist(permMeans(:,i),50);
end

% permuted distributions
plot(x,y,'color',ones(3,1)*.4)
xlabel('Mean difference (ms)')
ylabel('Count')

% observed difference
plot([1 1]*realdiff,get(gca,'ylim'),'r--','linew',4)

%% distribution of test statistics

figure(2), clf
subplot(211), hold on
hist(zval,30)
plot([1 1]*mean(zval),get(gca,'ylim'),'m:','linew',3)
xlabel('Z value'), ylabel('Count')
title([ 'Permuted Z values (average = ' num2str(mean(zval)) ')' ])


subplot(212), hold on
hist(pval,30)
plot([1 1]*mean(pval),get(gca,'ylim'),'m:','linew',3)
plot([1 1]*normcdf(-abs(mean(zval))),get(gca,'ylim'),'b:','linew',3)
xlabel('p value'), ylabel('Count')
title([ 'Permuted p values (average = ' num2str(normcdf(-abs(mean(zval)))) ')' ])

%% done.
