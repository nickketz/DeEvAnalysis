function outdata = deevGetEmerDep(cfg)
% function to calculate dependency using original Horner function and the
% emergent log files
%
% input:
%   cfg: config struct with optional fields
%        c: number of forced choice elements
%        filefiltstr: regexp filter on file listing, default is 'deev_36events_Sub[0-9]{1,2}\.txt'
%        badsubs: subject numbers to exclude from analysis
%        dir: struct containing the logs (def: logfiles)  
%        doplots: do plots or not (def 0)
%        ci: use confidence interval or not (def 1)
%
%
% output:
%   dep:    dependency with dimensions:
%           openXclosedObjXclosedAni X dataXindXdepXdepguess X locXper X cueXtarg  X subs
%   avgdep: dep averaged over dimensions 3 and 4, and rows 2 and 3,
%           i.e. canonical dependency matrix: openXclosed X dataXindXdepXdepguess X subs
%           
%       

%make MxN res for each subj
%
% how to structure N? I guess just do: 
% cues:
% 1:4 cue loc, 5:8 cue per, 9:12 cue obj, 13:16 cue ani
% targs:
% 1:4:16 targ loc, 2:4:16 targ per, 3:4:16 targ obj, 4:4:16 targ ani

%set defaults
if ~exist('cfg','var')          cfg = [];                end
if ~isfield(cfg,'filefiltstr')  cfg.filefiltstr = 'deev_36events_Sub[0-9]{1,2}\.txt';     end
if ~isfield(cfg,'badsubs')      cfg.badsubs = [];        end %which subjects to remove
if ~isfield(cfg,'dir')          cfg.dir = 'logs';        end %where are the log files?
if ~isfield(cfg,'c')            cfg.c = 6;               end
if ~isfield(cfg,'ci')           cfg.ci = 1;              end
if ~isfield(cfg,'doplots')      cfg.doplots = 0;         end



out = deevReadEmerLogs(cfg);
cfg_nt = cfg;
cfg_nt.filefiltstr = ['deev_nontarg_' cfg.filefiltstr(6:end)];
ntOut = deevReadEmerLogs(cfg_nt);

nsubs = size(out.res,3);
nEvts = size(out.res,1);
snums = regexp(out.lognames,'_Sub([0-9]+)\.txt$','tokens');
snums = cellfun(@(x) (str2num(x{1}{1})),snums);


%res = nan(nEvts,4*4,nsubs); % MxNxS = events X retrieval X subject
dep = nan(3,4,2,2,nsubs); % dependencies: openXclosed X dataXindXdepxdepgss X locXper X cueXtarg  X subs
E = nan(nEvts,2,2,nsubs); %episodic factor:  nEvts X locXper X cueXtarg  X subs 
condinds = cell(3,1);
condinds{1} = 1:2*nEvts/4;  
condinds{2} = 1 + 2*nEvts/4 : 3*nEvts/4; 
condinds{3} = 1 + 3*nEvts/4 : nEvts;
res = nan(size(out.res));

for isub = 1:nsubs
  
    tmpres = out.res(:,:,isub);
    %nan out nonsense retrievals
    %tmpres(condinds{1},[4 7 10 12 13 15]) = nan; %nan out loc-ani, per-obj, obj-per, obj-ani, ani-loc, ani-obj
    %tmpres(condinds{2}, [4 8 12 13 14 15]) = nan; %nan out loc-ani, per-ani, obj-ani, ani-loc, ani-per, ani-obj
    %tmpres(condinds{3}, [3 6 9 10 12 15]) = nan; %nan out loc-obj, per-obj, obj-loc, obj-per, obj-ani, ani-obj
    res(:,:,isub) = tmpres;

    %calc dependencies
    for icond = 1:3 %oepn, closedObj, closedAni
        
        switch icond
            case 1
                cuepair = [2 3; 5 8]; %[loc-per loc-obj; per-loc per-obj]
                targpair = [5 9; 2 14]; %[per-loc obj-loc; loc-per ani-per]
            case 2
                cuepair = [2 3; 5 7]; %[loc-per loc-obj; per-loc per-obj]
                targpair = [5 9; 2 10]; %[per-loc obj-loc; loc-per obj-per]
            case 3
                cuepair = [2 4; 5 8]; %[loc-per loc-ani; per-loc per-ani]
                targpair = [5 13; 2 14]; %[per-loc ani-loc; loc-per ani-per]
        end
        
        for iele = 1:2 %loc or per
            % cue
            [dep(icond,:,iele,1,isub), E(condinds{icond},iele,1,isub)] = deev_dependency(res(condinds{icond},:,isub),cuepair(iele,:),cfg.c);
            % targ
            [dep(icond,:,iele,2,isub), E(condinds{icond},iele,2,isub)] = deev_dependency(res(condinds{icond},:,isub),targpair(iele,:),cfg.c);            
        end
        
    end

end

avgdep = squeeze(mean(mean(dep,4),3));
avgdep = cat(1,avgdep(1,:,:), mean(avgdep(2:3,:,:)));

depdif = nan(size(avgdep,3),size(avgdep,2)-1,size(avgdep,1));
for icond = 1:size(avgdep,1)
    for idif = 2:size(avgdep,2)
        depdif(:,idif-1,icond) = avgdep(icond,1,:)-avgdep(icond,idif,:);
    end
end

oProb = squeeze(mean(res(1:nEvts/2,:,:),1));
oProb = reshape(oProb,4,[],nsubs);
oProb = permute(oProb,[2 1 3]);
altcProb = nanmean(res(1+nEvts/2:nEvts,:,:),1);
altcProb = reshape(altcProb,4,4,nsubs);
altcProb = permute(altcProb,[2 1 3]);

substr = strcat(repmat({'sub'},1,nsubs),strtrim(cellstr(num2str(snums'))'));

outdata.logdata = out;
outdata.subs = substr;

probDO = {{'cue-loc','cue-per','cue-obj','cue-ani'}, {'trg-loc','trg-per','trg-obj','trg-ani'}, substr};
outdata.oProb = oProb;
outdata.cProb = altcProb;
outdata.probDO = probDO;

outdata.dep = dep;
outdata.depDO = {{'openLoop','closedLoopObj','closedLoopAni'}, {'data','indp','dpnd','dpndg'}, {'loc','per'}, {'cue','targ'}, substr};

outdata.avgdep = avgdep;
outdata.avgdepDO = {{'openLoop','closedLoop'}, {'data','indp','dpnd','dpndg'}, substr};

outdata.avgdepdif = depdif;
outdata.avgdepdifDO = {substr,{'dIndp','dDep','dDepG'},{'openLoop','closedLoop'}};

outdata.res = res;
eles = {'loc','per','obj','ani'};
colstr = cell(1,size(tmpres,2));
for iele = 1:length(eles)
    for jele = 1:length(eles)
        colstr{(iele-1)*4 + jele} = [eles{iele} '-' eles{jele}];
    end
end
evtstr = strcat(repmat({'evt'},1,nEvts),strtrim(cellstr(num2str([1:nEvts]'))'));
outdata.resDO = {evtstr, colstr, substr};

outdata.resNT = ntOut.res;
outdata.resNTDO = outdata.resDO;

ntcl = ntOut.res(nEvts/2 + 1:end,:,:);
ntol = ntOut.res(1:nEvts/2,:,:);
outdata.ntacc = cat(1, nanmean(reshape(ntol,numel(ntol(:,:,1)),size(ntol,3))), nanmean(reshape(ntcl,numel(ntcl(:,:,1)),size(ntcl,3))));
outdata.ntaccDO = {{'openLoop','closedLoop'},substr};

%make ol and cl comparable a and b non-target vectors 
ab = squeeze(ntOut.res(:,2,:));
ac = cat(1,squeeze(ntOut.res(1:end-(nEvts/4),3,:)),squeeze(ntOut.res(3*nEvts/4 + 1:end,4,:)));
ba = squeeze(ntOut.res(:,5,:));
bc = cat(1,squeeze(ntOut.res(1:nEvts/2,8,:)),squeeze(ntOut.res(2*nEvts/4 + 1:3*nEvts/4,7,:)),squeeze(ntOut.res(3*nEvts/4 + 1:end,8,:)));
outdata.ntCueTarg = cat(3,ab,ac,ba,bc);
outdata.ntCueTargDO = {evtstr,substr, { 'loc-per','loc-oth','per-loc','per-oth'}};

outdata.E = E;
outdata.EDO = {evtstr,  {'loc','per'}, {'cue','targ'}, substr};


if cfg.doplots
    deevPlots(outdata,cfg);
end


        
        
        
    


