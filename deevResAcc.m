function acc = deevResAcc(in)
%
% function to calculate accuracy across cue and targ columns of deev 'res'
% matrix, seed deevGetLogDep or deevGetEmerDep for res matrix
%
% input:
%   in: res matrix with dimensions events x cue-targ array x subs
%
% ouput:
%   acc: 1 x subs accuracy collapsed across cue-targ array dimensions
%

acc = nanmean(reshape(in, numel(in(:,:,1)), size(in,3)));
