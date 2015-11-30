function outdata = deevGetEmerDep(cfg)
% function to calculate dependency using original Horner function and the
% emergent log files
%
% input:
%   cfg: config struct with optional fields
%        c: number of forced choice elements
%        filefiltstr: regexp filter on file listing, default is 'EntAssoc_*_Sub[0-9]{1,2}\.txt'
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
if ~isfield(cfg,'dir')          cfg.dir = 'logs';    end %where are the log files?
if ~isfield(cfg,'c')            cfg.c = 6;               end
if ~isfield(cfg,'ci')           cfg.ci = 1;              end
if ~isfield(cfg,'doplots')      cfg.doplots = 0;         end



out = deevReadEmerLogs(cfg);
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

outdata.E = E;
outdata.EDO = {evtstr,  {'loc','per'}, {'cue','targ'}, substr};


if cfg.doplots
    if cfg.ci
        crit = tinv(.975,nsubs-1);
    else
        crit = 1;
    end
    
    barstr = {'data','indp','dpnd','dpnd+g'};
    
    %plot dependency
    h = figure('color','white');
    mycolors = get(gca,'defaultAxesColorOrder');
    errorbar_groups(mean(avgdep,3)',crit*ste(avgdep,3)','bar_names',{'OpenLoop','ClosedLoop'},'bar_colors',mycolors,'FigID',h,...
        'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',5});
    legend(barstr,'location','best');
    title('Dependency frm log');
    ylim([min(min(mean(avgdep,3)))-.1 1]);
    set(gca,'fontsize',20);
    
    %plot by subj
    if size(avgdep,3)>1
    figure('color','white');
    hold on
    myconds = barstr;
%    mycolors = distinguishable_colors(length(myconds));
%    set(groot,'defaultAxesColorOrder',mycolors);
    lconds = {'openLoop','closedLoop'};
    mymin = 1; a = [];
    for ilcond = 1:2
        a(ilcond) = subplot(1,2,ilcond);
        tmp = squeeze(avgdep(ilcond,:,:));
        plot(tmp','.','markersize',30);
        hold on
        if ilcond == 1,        ylabel('Dependency','fontsize',18);       end
        xlabel('subject number','fontsize',18);
        tmp = tmp';
        for icond = 1:length(myconds)
            shadedErrorBar(0:size(tmp,1)+1,repmat(mean(tmp(:,icond)),[1 size(tmp,1)+2]),repmat(ste(tmp(:,icond)),[1 size(tmp,1)+2]),{'--','linewidth',2,'color',mycolors(icond,:),'markerfacecolor',mycolors(1,:)},1);
        end
        title(lconds{ilcond});
        set(gca,'fontsize',18);
        xlim([0 size(avgdep,3)]);
        box off
        mymin = min([mymin min(tmp)]);
    end
    legend(myconds,'fontsize',18,'location','southwest');
    for ilcond = 1:2,       ylim(a(ilcond),[mymin 1]);      end
    end
    
    
    %plot interaction
    h = figure('color','white');
    mycolors = get(gca,'defaultAxesColorOrder');
    errorbar_groups(squeeze(mean(depdif,1)),crit*squeeze(ste(depdif,1)),'bar_names',{'OpenLoop','ClosedLoop'},'bar_colors',mycolors(2:end,:),'FigID',h,...
        'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',5});
    legend({'dIndp','dDpnd','dDpnd+G'},'location','best');
    title('Delta Dependency');
    %ylim([min(depdif(:)) max(depdif(:))]);
    set(gca,'fontsize',20);
    
    
    %interaction by subject
    figure('color','white');
    hold on
    myconds = {'dInd','dDep','dDepG'};
    lconds = {'openLoop','closedLoop'};
    mymin = 1; mymax = 0; a = [];
    mycolors = mycolors(2:end,:);
    for ilcond = 1:2
        a(ilcond) = subplot(1,2,ilcond);
        tmp = squeeze(depdif(:,:,ilcond));
        set(gca,'ColorOrderIndex',3);
        for icond = 1:length(myconds)
            plot(tmp(:,icond),'.','markersize',30,'color',mycolors(icond,:));
            hold on
        end        
        if ilcond == 1,        ylabel('Delta Dependency','fontsize',18);       end
        xlabel('subject number','fontsize',18);
        %tmp = tmp';
        for icond = 1:length(myconds)
            shadedErrorBar(0:size(tmp,1)+1,repmat(mean(tmp(:,icond)),[1 size(tmp,1)+2]),repmat(crit*ste(tmp(:,icond)),...
                [1 size(tmp,1)+2]),{'--','linewidth',2,'color',mycolors(icond,:),'markerfacecolor',mycolors(1,:)},1);
        end
        title(lconds{ilcond});
        set(gca,'fontsize',18);
        xlim([0 size(avgdep,3)]);
        box off
        mymin = min([mymin min(tmp(:))]);
        mymax = max([mymax max(tmp(:))]);
    end
    legend(myconds,'fontsize',18,'location','southwest');
    for ilcond = 1:2,       ylim(a(ilcond),[mymin mymax]);      end
    
   
    %plot accuracy
    figure('color','white');
    mymin = min([min(min(nanmean(oProb,3))), min(min(nanmean(altcProb,3)))]);
    mymax = max([max(max(nanmean(oProb,3))), max(max(nanmean(altcProb,3)))]);
    h=subplot(1,2,1);
    imsc(nanmean(oProb,3), [mymin mymax], 'jet', [.5 .5 .5]);
    set(gca, 'DataAspectRatioMode', 'auto');
    set(h,'ytick',1:1:size(oProb,2));
    set(h,'yticklabel',probDO{1});
    set(h,'xtick', 1:size(oProb,1));
    set(h,'xticklabel',probDO{2});
    title('OpenLoop accuracy');
    set(gca,'fontsize',16);
    h=subplot(1,2,2);
    imsc(nanmean(altcProb,3), [mymin mymax], 'jet', [.5 .5 .5]); colorbar;
    set(gca, 'DataAspectRatioMode', 'auto');
    set(h,'ytick',1:1:size(altcProb,2));
    set(h,'yticklabel',probDO{1});
    set(h,'xtick', 1:size(altcProb,1));
    set(h,'xticklabel',probDO{2});
    title('ClosedLoop accuracy');
    set(gca,'fontsize',16);
    
    
    
    
    %plot by type
    h = figure('color','white');
    mycolors = distinguishable_colors(8,{'w','k'});
    hold on
    mymin = 1;
    a = [];
    stops = [-1.5 -.5 .5 1.5];
    for iaxis = 1:2
        a(iaxis) = subplot(1,2,iaxis);
        if iaxis == 1, idim = 2; else idim = 1; end %hack to fix mean dimenison mismatch
        omu = squeeze(nanmean(oProb,idim)); if size(omu,1)~= 4, omu = omu'; end
        cmu = squeeze(nanmean(altcProb,idim)); if size(cmu,1)~=4, cmu = cmu'; end
        mu = cat(3,omu,cmu);
        mu = permute(mu, [3 1 2]);
        
        [xtick,hb,he]=errorbar_groups(mean(mu,3)', crit*ste(mu,3)','bar_names',{'OpenLoop','ClosedLoop'},'FigID', h, 'AxID', a(iaxis),...
            'bar_colors',mycolors(1+((iaxis-1)*size(mycolors,1)/2):(iaxis*size(mycolors,1)/2),:),...
            'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',5});
        
        legend(probDO{iaxis});
        hold on
        vals = mean(mu,3);
        for ib = 1:length(hb)
            plot(repmat(xtick(1)+stops(ib),size(mu,3),2),squeeze(mu(1,ib,:)),'k.','markersize',20)
            plot(repmat(xtick(2)+stops(ib),size(mu,3),2),squeeze(mu(2,ib,:)),'k.','markersize',20)
        end
            
        
        if iaxis == 1,          ylabel('Accuracy');        end
        set(a(iaxis),'fontsize',20);
        if size(mu,3)>1
            mymin = min([mymin min(mean(mu,3)-(crit*ste(mu,3)))]);
        else
            mymin = min([mymin min(mu(:)-.1)]);
        end
    end
    for iaxis = 1:2
        ylim(a(iaxis),[mymin 1]);
    end
    
end


        
        
        
    


