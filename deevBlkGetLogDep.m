function blkdata = deevBlkGetLogDep(cfg)
% do Dependency analysis by block
% for extended deev experiment with multiple blocks, analysis does
% dependency analysis for each block independently using the original
% deevGetLogDep function
%
% input
%   cfg: config struct with optional fields
%        c: number of forced choice elements
%        filefiltstr: regexp filter on file listing, default is 'deevEEGExp_0EEG_4blocks_36events_Sub[0-9]{1,2}\.txt'
%        badsubs: subject numbers to exclude from analysis
%        dir: struct containing the logs (def: logfiles)
%        doplots: do plots or not (def 0)
%        ci: use confidence interval or not (def 1)
%
% output:
%   blkdata: output struct with fields
%       dep: dependency for each block with dimensions labeld in depDO:
%               openXclosedObjXclosedAni X dataXindXdepXdepguess X locXper X
%               cueXtarg  X subs X blk
%       avgdep: dep averaged over dimensions 3 and 4, and rows 2 and 3,
%               i.e. canonical dependency matrix: openXclosed X
%               dataXindXdepXdepguess X subs; dimension labeled in avgdepDO
%       oProb: open loop accuracy by block, dimesions labeled in oProbDO
%       cProb: closed loop accuracy by block, dimensions in cProbDO
%       res: correct response binary matrix by block, dimensions in resDO
%       E: episodic factor by block, dimension labls in EDO
%       
%       
%
%

%set defaults
if ~exist('cfg','var')          cfg = [];                end
if ~isfield(cfg,'filefiltstr')  cfg.filefiltstr = 'deevEEGExp_0EEG_4blocks_36events_Sub[0-9]{1,2}\.txt';     end
if ~isfield(cfg,'badsubs')      cfg.badsubs = [];        end %which subjects to remove
if ~isfield(cfg,'dir')          cfg.dir = 'logfiles';    end %where are the log files?
if ~isfield(cfg,'c')            cfg.c = 6;               end
if ~isfield(cfg,'ci')           cfg.ci = 1;              end
if ~isfield(cfg,'doplots')      cfg.doplots = 0;         end



out = deevReadlogs(cfg);
nEvts = max(out.data.Trial(:));
nBlks = max(out.data.ExpBlock(:));

blks = unique(out.data.ExpBlock(:,1),'stable');
blkcnt = 0;
for iblk = blks'
    blkcnt = blkcnt + 1;
    %get inds
    blkind = out.data.ExpBlock(:,1) == iblk;
    if sum(blkind) ~= 3*nEvts*3, error('blkind doesn''t macth nEvts*3'); end
    p2ind = (out.data.Phase(:,1) == 2) & blkind;
    if sum(p2ind) ~= 2*sum(blkind)/3, error('p2ind doesn''t 2/3 of total block trials'); end
    
    %get vars
    crrct = out.data.Correct(p2ind,:);
    ord = out.data.Order(p2ind,:);
    stmtp = out.data.StimType(p2ind,:);
    evt = out.data.Trial(p2ind,:);    
    nsubs = size(crrct,2);
    
    res = nan(nEvts,4*4,nsubs); % MxNxS = events X retrieval X subject
    dep = nan(3,4,2,2,nsubs); % dependencies: openXclosedObjXclosedAni X dataXindXdepxdepgss X locXper X cueXtarg  X subs
    E = nan(nEvts,2,2,nsubs); %episodic factor:  nEvts X locXper X cueXtarg  X subs
    condinds = cell(3,1);
    condinds{1} = 1:2*nEvts/4;
    condinds{2} = 1 + 2*nEvts/4 : 3*nEvts/4;
    condinds{3} = 1 + 3*nEvts/4 : nEvts;
    
    
    for isub = 1:nsubs
        
        %get response accuracy
        for itrl = 1:size(crrct,1)
            %get event number
            tmpEvt = evt(itrl,isub);
            
            %sanity check on logs
            if tmpEvt > nEvts/2
                if ~strcmp(stmtp,'CL')
                    error('closed loop stim types are not in the second half of events\nsub: %i, row: %i, trial: %i',isub,itrl,tmpEvt);
                end
            else
                if ~strcmp(stmtp,'OL')
                    error('open loop stim types are not in the first half of events\nsub: %i, row: %i, trial: %i',isub,itrl,tmpEvt);
                end
            end
            
            %get cue and targ stim categories
            tmpOrd = regexp(ord(itrl,isub),'-','split');
            tmpOrd = cellfun(@str2num,tmpOrd{1});
            %replace obj with ani for last quater of closed loop trials
            %don't need this anymore, fixed in experiment script
            %if tmpEvt > 3*nEvts/4,      tmpOrd(tmpOrd==3) = 4;          end
            
            %calc column index for this trial
            colInd = ((tmpOrd(1)-1)*4)+1 + (tmpOrd(2)-1);
            
            %assign accuracy to res matrix
            res(tmpEvt,colInd,isub) = crrct(itrl,isub);
            
        end
        
        
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
    cProb = nanmean(res(1+nEvts/2:nEvts,:,:),1);
    cProb = reshape(cProb,4,4,nsubs);
    cProb = permute(cProb,[2 1 3]);
    
    %store in block vars
    if blkcnt == 1
        blkoProb = oProb;
        blkcProb = cProb;
        blkAvgDep = avgdep;
        blkDep = dep;
        blkE = E;
        blkRes = res;
        blkDepDif = depdif;
    else
        blkoProb = cat(length(size(oProb))+1,blkoProb,oProb);
        blkcProb = cat(length(size(cProb))+1,blkcProb,cProb);
        blkAvgDep = cat(length(size(avgdep))+1,blkAvgDep,avgdep);
        blkDep = cat(length(size(dep))+1,blkDep,dep);
        blkE = cat(length(size(E))+1,blkE,E);
        blkRes = cat(length(size(res))+1,blkRes,res);
        blkDepDif = cat(length(size(depdif))+1,blkDepDif,depdif);
    end
        
end

%construct dimension labels
substr = strcat(repmat({'sub'},1,nsubs),strtrim(cellstr(num2str(out.data.Subj(1,:)'))'));
blkstr = strcat(repmat({'blk'},1,nBlks),strtrim(cellstr(num2str([1:nBlks]'))'));
evtstr = strcat(repmat({'evt'},1,nEvts),strtrim(cellstr(num2str([1:nEvts]'))'));
eles = {'loc','per','obj','ani'};
colstr = cell(1,size(res,2));
for iele = 1:length(eles)
    for jele = 1:length(eles)
        colstr{(iele-1)*4 + jele} = [eles{iele} '-' eles{jele}];
    end
end

blkdata.logdata = out;
blkdata.subs = substr;

blkdata.oProb = blkoProb;
blkdata.cProb = blkcProb;
blkdata.probDO = {{'cue-loc','cue-per','cue-obj','cue-ani'}, {'trg-loc','trg-per','trg-obj','trg-ani'}, substr, blkstr};

blkdata.dep = blkDep;
blkdata.depDO = {{'openLoop','closedLoopObj','closedLoopAni'}, {'data','indp','dpnd','dpndg'}, {'loc','per'}, {'cue','targ'}, substr, blkstr};

blkdata.avgdep = blkAvgDep;
blkdata.avgdepDO = {{'openLoop','closedLoop'}, {'data','indp','dpnd','dpndg'}, substr, blkstr};

blkdata.avgdepdif = blkDepDif;
blkdata.avgdepdifDO = {substr,{'dIndp','dDep','dDepG'},{'openLoop','closedLoop'},blkstr};

blkdata.res = blkRes;
blkdata.resDO = {evtstr, colstr, substr, blkstr};

blkdata.E = blkE;
blkdata.EDO = {evtstr, {'data','indp','dpnd','dpndg'}, {'loc','per'}, {'cue','targ'}, substr, blkstr};

%create block averages for plotting
avgdep = squeeze(mean(blkdata.avgdep,4));
oProb = squeeze(mean(blkdata.oProb,4));
cProb = squeeze(mean(blkdata.cProb,4));
probDO = blkdata.probDO{1:end-1};
if cfg.doplots
    
    deevPlots(blkdata,cfg)
%     if cfg.ci
%         crit = tinv(.975,nsubs-1);
%     else
%         crit = 1;
%     end
%     
%     barstr = {'data','indp','dpnd','dpnd+g'};
%     
%     %plot dependency
%     h = figure('color','white');
%     mycolors = get(gca,'defaultAxesColorOrder');
%     errorbar_groups(mean(avgdep,3)',crit*ste(avgdep,3)','bar_names',{'OpenLoop','ClosedLoop'},'bar_colors',mycolors,'FigID',h,...
%         'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',5});
%     legend(barstr,'location','best');
%     title('Dependency frm log');
%     ylim([min(min(mean(avgdep,3)))-.1 1]);
%     set(gca,'fontsize',20);
%     
%     %plot by subj
%     figure('color','white');
%     hold on
%     myconds = barstr;
%     lconds = {'openLoop','closedLoop'};
%     mymin = 1; a = [];
%     for ilcond = 1:2
%         a(ilcond) = subplot(1,2,ilcond);
%         tmp = squeeze(avgdep(ilcond,:,:));
%         plot(tmp','.','markersize',30);
%         hold on
%         if ilcond == 1,        ylabel('Dependency','fontsize',18);       end
%         xlabel('subject number','fontsize',18);
%         tmp = tmp';
%         for icond = 1:length(myconds)
%             shadedErrorBar(0:size(tmp,1)+1,repmat(mean(tmp(:,icond)),[1 size(tmp,1)+2]),repmat(ste(tmp(:,icond)),[1 size(tmp,1)+2]),{'--','linewidth',2,'color',mycolors(icond,:),'markerfacecolor',mycolors(1,:)},1);
%         end
%         title(lconds{ilcond});
%         set(gca,'fontsize',18);
%         xlim([0 size(avgdep,3)]);
%         box off
%         mymin = min([mymin min(tmp)]);
%     end
%     legend(myconds,'fontsize',18,'location','southwest');
%     for ilcond = 1:2,       ylim(a(ilcond),[mymin 1]);      end
%     
%     
%     
%     
%     %plot accuracy
%     figure('color','white');
%     mymin = min([min(min(nanmean(oProb,3))), min(min(nanmean(cProb,3)))]);
%     mymax = max([max(max(nanmean(oProb,3))), max(max(nanmean(cProb,3)))]);
%     h=subplot(1,2,1);
%     imsc(nanmean(oProb,3), [mymin mymax], 'jet', [.5 .5 .5]);
%     set(gca, 'DataAspectRatioMode', 'auto');
%     set(h,'ytick',1:1:size(oProb,2));
%     set(h,'yticklabel',probDO{1});
%     set(h,'xtick', 1:size(oProb,1));
%     set(h,'xticklabel',probDO{2});
%     title('OpenLoop accuracy');
%     set(gca,'fontsize',16);
%     h=subplot(1,2,2);
%     imsc(nanmean(cProb,3), [mymin mymax], 'jet', [.5 .5 .5]); colorbar;
%     set(gca, 'DataAspectRatioMode', 'auto');
%     set(h,'ytick',1:1:size(cProb,2));
%     set(h,'yticklabel',probDO{1});
%     set(h,'xtick', 1:size(cProb,1));
%     set(h,'xticklabel',probDO{2});
%     title('ClosedLoop accuracy');
%     set(gca,'fontsize',16);
%     
%     %plot by type
%     h = figure('color','white');
%     mycolors = distinguishable_colors(8,{'w','k'});
%     hold on
%     mymin = 1;
%     a = [];
%     for iaxis = 1:2
%         a(iaxis) = subplot(1,2,iaxis);
%         if iaxis == 1, idim = 2; else idim = 1; end %hack to fix mean dimenison mismatch
%         omu = squeeze(nanmean(oProb,idim));
%         cmu = squeeze(nanmean(cProb,idim));
%         mu = cat(3,omu,cmu);
%         mu = permute(mu, [3 1 2]);
%         errorbar_groups(mean(mu,3)', crit*ste(mu,3)','bar_names',{'OpenLoop','ClosedLoop'},'FigID', h, 'AxID', a(iaxis),...
%             'bar_colors',mycolors(1+((iaxis-1)*size(mycolors,1)/2):(iaxis*size(mycolors,1)/2),:),...
%             'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',5});
%         legend(probDO{iaxis});
%         if iaxis == 1,          ylabel('Accuracy');        end
%         set(a(iaxis),'fontsize',20);
%         mymin = min([mymin min(mean(mu,3)-(crit*ste(mu,3)))]);
%     end
%     for iaxis = 1:2
%         ylim(a(iaxis),[mymin 1]);
%     end
%     
end


