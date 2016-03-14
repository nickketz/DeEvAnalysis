function [avgtestsim,avgtrainsim,twin] = deevCrosValW(in,slide,data,nfold)
% perform cross validation on selected time/freq window
% returns max dot product for each trial in each train and test fold
%
% input:
%   in: 4 element vector: tstart, tend, fstart, fend
%   slide: temporal slide to evaluate w on
%   data: ft pow struct data to derive W from, must have subselected chans
%   to match W
%   nfold: number of cross validation folds to run, -1 = leave one out
%
% output:
%   avgtrainsim: nfold x nslide; normpord(w,train_data) slides for a
%   given fold, then averaged across trials
%   test: ntrial x nslide; normprod(w,test_data) slide across time for each
%   trial
%   twin: nslide x 2; start and end of time windows used in slides
%


%break out input vector
s = in(1); tdur = in(2); l = in(3); fdur = in(4);


%grab data to make template
[~,~,s] = findnearest(s,data.time);
[~,~,l] = findnearest(l,data.freq);
%tdur = ((max(data.time)-s)*tpct);
%fdur = ((max(data.freq)-l)*fpct);

tsel = data.time>=s & data.time<=s+tdur;
fsel = data.freq>=l & data.freq<=l+fdur;

%make hold-out set
ntrials = size(data.powspctrm,1);
if nfold == -1 || nfold > ntrials
    nval = ntrials;
else
    nval = nfold;
end

pctOut = 1/nval;
if nfold > 0
    nOut = floor(pctOut*ntrials);
else
    nOut = 1;
end
testinds = zeros(nval,ntrials);
for ival = 1:nval, testinds(ival,(ival-1)*nOut + 1:ival*nOut)=1; end
testinds = testinds(:,randperm(ntrials)); testinds = logical(testinds(:,1:ntrials));
traininds = ~testinds;


%compare tempalte with observations
slide = slide+mod(slide,diff(data.time(1:2))); %make slide divisible by time grains
nsamp = floor((diff(data.time(1:2))*length(data.time)-tdur)/slide);
testsim = nan(nOut,nsamp);
avgtestsim = [];
avgtrainsim = [];
fsamp = fsel;
twin = nan(nsamp,2);

fprintf('\nrunning %d folds: 00%%',nval);
for ival = 1:nval
    fprintf('\b\b\b%02d%%',round(100*ival/nval));
    trainsel = traininds(ival,:);
    testsel = testinds(ival,:);
    trainsim = nan(sum(trainsel),nsamp);
    
    w = squeeze(mean(data.powspctrm(trainsel,:,fsel,tsel),1));
    w = reshape(w,1,numel(w(:)));
    
    for isamp = 0:nsamp
        [~,i] = findnearest(isamp*slide + data.time(1),data.time);
        tsamp = i:i+sum(tsel)-1;
        twin(isamp+1,:) = [data.time(i), data.time(tsamp(end))];
        
        trainx = data.powspctrm(trainsel,:,fsamp,tsamp);
        trainx = reshape(trainx,size(trainx,1),numel(trainx(1,:,:,:)));
        
        testx = data.powspctrm(testsel,:,fsamp,tsamp);
        testx = reshape(testx,size(testx,1),numel(testx(1,:,:,:)));
        
        
        if size(testx)~=size(w),         error('mismatch sample size');      end
        
        testsim(:,isamp+1) = normprod(w,testx')';
        trainsim(:,isamp+1) = normprod(w,trainx')';
    end
%    avgtestsim = cat(2,avgtestsim, max(testsim,[],2));
%    avgtrainsim = cat(2,avgtrainsim, mean(max(trainsim,[],2)));
    avgtestsim = cat(1,avgtestsim, testsim);
    avgtrainsim = cat(1,avgtrainsim, mean(trainsim,1));
end
fprintf('\b\b\b\bdone');
