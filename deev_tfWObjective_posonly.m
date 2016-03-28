function score = deev_tfWObjective(in)
% objective function to evaluate time frequency template in optimatzation
% search
% template W is the average overtraining trials in the timefreq window
% specified by 'in'
%
% this template is evaluated in a 4 fold cross validation based on the
% normprod of W with the hold-out trials in a sliding window manner
% starting with t=0 of the hold-out trials and evaluating at at t+slide
% until W no longer fits
%
% score is returned as the inverse max normprod across slides for each
% trial, then averaged across trials
% (i.e. max of slides, averaged over trials)
%
%
% input:
%   in: 4 parameter vector:
%       x1=start time, x2=duration time window after start to use
%       x3=start freq, x4=length of freq window after start to use
%
%       note: kinda clunky pct parameters used to avoid complications of
%       end time/freq can't be less then start time/freq
%
% output:
%   score: evaluation score is one over the average over trials of the max
%   normprod over slides.

%grab variables from base workspace
dpos = evalin('base','dpos');
%dneg = evalin('base','dneg');
chaninds = evalin('base','chaninds');
slide = evalin('base','slide');

%break out input vector
s = in(1); tdur = in(2); l = in(3); fdur = in(4);

%check if params make sense, return Inf is not
if s+tdur>dpos.time(end) 
    score = Inf;
    return
end

if l+fdur>dpos.freq(end)
    score = Inf;
    return
end



%grab data to make template
[~,~,s] = findnearest(s,dpos.time);
[~,~,l] = findnearest(l,dpos.freq);
%tdur = ((max(dpos.time)-s)*tpct);
%fdur = ((max(dpos.freq)-l)*fpct);

tsel = dpos.time>=s & dpos.time<=s+tdur;
fsel = dpos.freq>=l & dpos.freq<=l+fdur;

%make hold-out set
nval = 4;
ntrials = size(dpos.powspctrm,1);
pctOut = 1/nval;
nOut = floor(pctOut*ntrials);
testinds = zeros(nval,ntrials);
for ival = 1:nval, testinds(ival,(ival-1)*nOut + 1:ival*nOut)=1; end
testinds = testinds(:,randperm(ntrials)); testinds = logical(testinds(:,1:ntrials));
traininds = ~testinds;


%compare tempalte with observations
nsamp = floor(max(dpos.time-tdur)/slide);
possim = nan(nOut,nsamp);
avgpossim = [];
%negsim = nan(size(dneg.powspctrm,1),nsamp);
%avgnegsim = [];
fsamp = fsel;

for ival = 1:nval
    
    trainsel = traininds(ival,:);
    testsel = testinds(ival,:);
    
    w = squeeze(mean(dpos.powspctrm(trainsel,chaninds,fsel,tsel),1));
    w = reshape(w,1,numel(w(:)));
    
    for isamp = 0:nsamp
        [~,i] = findnearest(isamp*slide,dpos.time);
        tsamp = i:i+sum(tsel)-1;
        
        posx = dpos.powspctrm(testsel,chaninds,fsamp,tsamp);
        posx = reshape(posx,size(posx,1),numel(posx(1,:,:,:)));
        
%         negx = dneg.powspctrm(:,chaninds,fsamp,tsamp);
%         negx = reshape(negx,size(negx,1),numel(negx(1,:,:,:)));
        
        if size(posx)~=size(w),         error('mismatch sample size');      end
        
        possim(:,isamp+1) = normprod(w,posx')';
%        negsim(:,isamp+1) = normprod(w,negx')';
    end
    avgpossim = cat(2,avgpossim, max(possim,[],2));
%    avgnegsim = cat(2,avgnegsim, max(negsim,[],2));
%     if ival == 1
%         avgpossim = max(possim);
%         avgnegsim = negsim;
%     else
%         avgpossim = avgpossim + max(possim);
%         avgnegsim = avgnegsim + negsim;
%     end
    
end

%avgpossim = avgpossim./nval;
%avgnegsim = avgnegsim./nval;

% define scale param to make min focus first on get high avgpossim score,
% then getting avgnegsim as close to 0 as possible
%mean(avgpossim) = .5 => exp portion = .1, i.e. balanced influence between
%neg and pos cases 
scale = 9; 
score = exp(-scale*mean(avgpossim(:)));
%score = exp(-scale*mean(avgpossim(:))) + abs(mean(avgnegsim(:)));
%score = exp(-scale*mean(max(avgpossim))) + abs(mean(max(avgnegsim)));
%score = 1/(mean(max(avgpossim))) + (1-abs(mean(max(avgnegsim))));


