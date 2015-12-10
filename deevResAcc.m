function [acc, oste] = deevResAcc(in)
%
% function to calculate accuracy across cue and targ columns of deev 'res'
% matrix, seed deevGetLogDep or deevGetEmerDep for res matrix
%
% input:
%   in: res matrix with dimensions events x cue-targ array x subs
%
% ouput:
%   acc: blocks x subs accuracy collapsed across cue-targ array dimensions
%

acc = [];
oste = [];
if size(in,4)==4
    %do by block
    for i = 1:size(in,4)
        tmp = in(:,:,:,i);
        acc = cat(1,acc, nanmean(reshape(tmp, numel(tmp(:,:,1)), size(tmp,3))));
        oste = cat(1,oste, nanste(reshape(tmp, numel(tmp(:,:,1)), size(tmp,3))));
    end
else
    acc = nanmean(reshape(in, numel(in(:,:,1)), size(in,3)));
    oste = nanste(reshape(in, numel(in(:,:,1)), size(in,3)));
end