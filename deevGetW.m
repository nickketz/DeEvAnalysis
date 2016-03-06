function [w,t,f,chans] = deevGetW(x)
% create tempalte W from categorical stim search
% dpos and chaninds still need to be in base workspace from search  
% 
% input:
%   x = optimal 4 parameter vector from search
%
% output
%   w = optimal template in vector form
%   t = optimal time points used to make W
%   f = optimal freq points use to make W
%   chans = channels used to make W

dpos = evalin('base','dpos');
chaninds = evalin('base','chaninds');

s = x(1); tdur = x(2); l = x(3); fdur = x(4);

%grab data to make template
[~,~,s] = findnearest(s,dpos.time);
[~,~,l] = findnearest(l,dpos.freq);
%tdur = ((max(dpos.time)-s)*tpct);
%fdur = ((max(dpos.freq)-l)*fpct);

tsel = dpos.time>=s & dpos.time<=s+tdur;
fsel = dpos.freq>=l & dpos.freq<=l+fdur;

t = dpos.time(tsel);
f = dpos.freq(fsel);
chans = dpos.label(chaninds);

w = squeeze(mean(dpos.powspctrm(:,chaninds,fsel,tsel),1));
w = reshape(w,1,numel(w(:)));
