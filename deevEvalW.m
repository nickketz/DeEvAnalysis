function [sim, tsamps] = deevEvalW(w,t,f,chaninds,data,slide)
%
% function to evaluate template w on a give condition data
% 
% input;
%   w: vector, optimal template to use in evaluations
%   t: vector, optimal time points used
%   f: vector, optimal freq points used
%   chaninds: logical, subsample of channels used in W
%   data: trial level ft pow struct to evalute W on
%   slide: temporal distance between evaluations in secs (def: 0.1s)
%
% output:
%   sim: evals x trials similarity score for each trial, and each 'slide' in data
%   tsamps: evals x twindow time samples used in each eval of sim
%

if ~exist('slide','var');           slide = .1;                         end

s = min(t); e = max(t);
l = min(f); u = max(f);

%grab data to make template
[~,~,s] = findnearest(s,data.time);
[~,~,l] = findnearest(l,data.freq);
tdur = t(end)-t(1);
fdur = f(end)-f(1);

tsel = data.time>=s & data.time<=s+tdur;
fsel = data.freq>=l & data.freq<=l+fdur;

if sum(data.time(tsel)~=t)>0 || sum(data.freq(fsel)~=f)>0
    error('selected time/freq dimensions don''t match optimal');
end

%compare tempalte with observations
nsamp = floor(max(data.time-tdur)/slide);
ntrials = size(data.trialinfo,1);
sim = nan(ntrials,nsamp);
fsamp = fsel;
tsamps = nan(ntrials,length(t));

for isamp = 0:nsamp
    [~,i] = findnearest(isamp*slide,data.time);
    tsamp = i:i+sum(tsel)-1;
    tsamps(isamp+1,:) = data.time(tsamp);
    
    x = data.powspctrm(:,chaninds,fsamp,tsamp);
    x = reshape(x,size(x,1),numel(x(1,:,:,:)));
    
    if size(x)~=size(w),         error('mismatch sample size');      end
    
    sim(:,isamp+1) = normprod(w,x')';
end
    
