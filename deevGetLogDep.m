function outdata = deevGetLogDep(cfg)
% function to calculate dependency using original Horner function and the
% log files as opposed to the mat files
%
% input:
%   cfg: config struct with optional fields
%        c: number of forced choice elements
%        filefiltstr: regexp filter on file listing, default is 'deevEEGExp_0EEG_4blocks_36events_Sub[0-9]{1,2}\.txt'
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
if ~isfield(cfg,'filefiltstr')  cfg.filefiltstr = 'deevEEGExp_0EEG_4blocks_36events_Sub[0-9]{1,2}\.txt';     end
if ~isfield(cfg,'badsubs')      cfg.badsubs = [];        end %which subjects to remove
if ~isfield(cfg,'dir')          cfg.dir = 'logfiles';    end %where are the log files?
if ~isfield(cfg,'c')            cfg.c = 6;               end
if ~isfield(cfg,'ci')           cfg.ci = 1;              end
if ~isfield(cfg,'doplots')      cfg.doplots = 0;         end



out = deevReadlogs(cfg);
nEvts = max(out.data.Trial(:));

p2ind = out.data.Phase(:,1) == 2;
crrct = out.data.Correct(p2ind,:);
ord = out.data.Order(p2ind,:);
stmtp = out.data.StimType(p2ind,:);
evt = out.data.Trial(p2ind,:);

nsubs = size(crrct,2);
nrets = size(crrct,1);

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
        if tmpEvt > 3*nEvts/4,      tmpOrd(tmpOrd==3) = 4;          end 
        
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

oProb = squeeze(mean(res(1:nEvts/2,:,:),1));
oProb = reshape(oProb,4,[],nsubs);
oProb = permute(oProb,[2 1 3]);
altcProb = nanmean(res(1+nEvts/2:nEvts,:,:),1);
altcProb = reshape(altcProb,4,4,nsubs);
altcProb = permute(altcProb,[2 1 3]);

depdif = nan(size(avgdep,3),size(avgdep,2)-1,size(avgdep,1));
for icond = 1:size(avgdep,1)
    for idif = 2:size(avgdep,2)
        depdif(:,idif-1,icond) = avgdep(icond,1,:)-avgdep(icond,idif,:);
    end
end

substr = strcat(repmat({'sub'},1,nsubs),strtrim(cellstr(num2str(out.data.Subj(1,:)'))'));

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
colstr = cell(1,size(res,2));
for iele = 1:length(eles)
    for jele = 1:length(eles)
        colstr{(iele-1)*4 + jele} = [eles{iele} '-' eles{jele}];
    end
end
evtstr = strcat(repmat({'evt'},1,nEvts),strtrim(cellstr(num2str([1:nEvts]'))'));
outdata.resDO = {evtstr, colstr, substr};

outdata.E = E;
outdata.EDO = {evtstr, {'data','indp','dpnd','dpndg'}, {'loc','per'}, {'cue','targ'}, substr};


if cfg.doplots
    
    deevPlots(outdata,cfg);
%     
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
% %    mycolors = distinguishable_colors(length(myconds));
% %    set(groot,'defaultAxesColorOrder',mycolors);
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
%     mymin = min([min(min(nanmean(oProb,3))), min(min(nanmean(altcProb,3)))]);
%     mymax = max([max(max(nanmean(oProb,3))), max(max(nanmean(altcProb,3)))]);
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
%     imsc(nanmean(altcProb,3), [mymin mymax], 'jet', [.5 .5 .5]); colorbar;
%     set(gca, 'DataAspectRatioMode', 'auto');
%     set(h,'ytick',1:1:size(altcProb,2));
%     set(h,'yticklabel',probDO{1});
%     set(h,'xtick', 1:size(altcProb,1));
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
%         cmu = squeeze(nanmean(altcProb,idim));
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


        
        
        
    


