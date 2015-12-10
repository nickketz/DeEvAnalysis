function deevPlots(in,cfg)
%
%function to produce behavioral plots for deev experiments
%
% input:
%   in: data struct from deevGetLogDep/deevBlkGetLogDep
%   cfg: config struct
%       ci = bool confidence interval? def=1
%       blk = array of which blocks to do plots on, or 'avg' implies an
%             average across blocks, def='avg'
%       mdl = cell array of which models to plot, default {'ind','dpn','dpng'}
%       blkplts = bool, do block plots of accuracy, dependency and ol vs cl
%       delta dep, def=0
%
%
%

%set defaults
if ~exist('cfg','var')      cfg = [];                                   end
if ~isfield(cfg,'ci')       cfg.ci = 1;                                 end
if ~isfield(cfg,'blk')      cfg.blk = 'avg';                            end
if ~isfield(cfg,'mdls')     cfg.mdls = {'data','indp','dpnd','dpndg'};  end
if ~isfield(cfg,'blkplts')  cfg.blkplts = 0;                            end


%is in from the blk version of the study?
blkstr = '';
if length(size(in.dep))==6
    
    if strcmp(cfg.blk, 'avg')
        %create block averages for plotting
        avgdep = squeeze(mean(in.avgdep,4));
        avgdepdif = squeeze(mean(in.avgdepdif,4));
        oProb = squeeze(mean(in.oProb,4));
        cProb = squeeze(mean(in.cProb,4));
        probDO = in.probDO(1:end-1);
        blkstr = 'avg';
    else
        %sub sample specific block
        avgdep = in.avgdep(:,:,:,cfg.blk);
        avgdepdif = in.avgdepdif(:,:,:,cfg.blk);
        oProb = in.oProb(:,:,:,cfg.blk);
        cProb = in.cProb(:,:,:,cfg.blk);
        probDO = in.probDO(1:end-1);
        blkstr = num2str(cfg.blk);
    end
    
else %not in blocks
    avgdep = in.avgdep;
    avgdepdif = in.avgdepdif;
    oProb = in.oProb;
    cProb = in.cProb;
    probDO = in.probDO;
end
    
nsubs = length(in.logdata.lognames);

mdls = {'data','indp','dpnd','dpndg'};
mdlinds = find(ismember(cfg.mdls,mdls));

if cfg.ci
    crit = tinv(.975,nsubs-1);
else
    crit = 1;
end

barstr = {'data','indp','dpnd','dpnd+g'};


%plot dependency
h = figure('color','white','name',blkstr);
mycolors = get(gca,'defaultAxesColorOrder');
tmp = avgdep(:, mdlinds, :);
errorbar_groups(mean(tmp,3)',crit*ste(tmp,3)','bar_names',{'OpenLoop','ClosedLoop'},'bar_colors',mycolors,'FigID',h,...
    'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',5});
legend(barstr(mdlinds),'location','best');
title('Dependency');
%ylim([min(min(mean(tmp,3)))-.1 1]);
ylim([.6 1]);
set(gca,'fontsize',20);

%plot by subj
if size(avgdep,3)>1
    figure('color','white','name',blkstr);
    hold on
    myconds = barstr;
    %    mycolors = distinguishable_colors(length(myconds));
    %    set(groot,'defaultAxesColorOrder',mycolors);
    lconds = {'openLoop','closedLoop'};
    mymin = 1; a = [];
    for ilcond = 1:2
        a(ilcond) = subplot(1,2,ilcond);
        tmp = squeeze(avgdep(ilcond,mdlinds,:));
        plot(tmp','.','markersize',30);
        hold on
        if ilcond == 1,        ylabel('Dependency','fontsize',18);       end
        xlabel('subject number','fontsize',18);
        tmp = tmp';
        for icond = 1:length(myconds(mdlinds))
            shadedErrorBar(0:size(tmp,1)+1,repmat(mean(tmp(:,icond)),[1 size(tmp,1)+2]),repmat(ste(tmp(:,icond)),[1 size(tmp,1)+2]),{'--','linewidth',2,'color',mycolors(icond,:),'markerfacecolor',mycolors(1,:)},1);
        end
        title(lconds{ilcond});
        set(gca,'fontsize',18);
        xlim([0 size(avgdep,3)]);
        box off
        mymin = min([mymin min(tmp)]);
    end
    legend(myconds(mdlinds),'fontsize',18,'location','southwest');
    mymin = min(avgdep(:)); mymax = max(avgdep(:));
    for ilcond = 1:2,       ylim(a(ilcond),[mymin mymax]);      end
end


%plot interaction
h = figure('color','white','name',blkstr);
mycolors = get(gca,'defaultAxesColorOrder');
[bar,hb,he] = errorbar_groups(squeeze(mean(avgdepdif,1)),crit*squeeze(ste(avgdepdif,1)),'bar_names',{'OpenLoop','ClosedLoop'},'bar_colors',mycolors(2:end,:),'FigID',h,...
    'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',5});
legend({'dIndp','dDpnd','dDpnd+G'},'location','best');
title('Delta Dependency');
%ylim([min(depdif(:)) max(depdif(:))]);
set(gca,'fontsize',20);
%plot significance?
mu = squeeze(mean(avgdepdif,1));
for ibar = 1:length(bar)
    [sig,p] = ttest(avgdepdif(:,:,ibar));
    for isig = 1:length(sig)
        if sig(isig)
            dif = [-1 0 1];
            ptext = '*';
            if p(isig)<0.01     ptext = [ptext '*'];    end
            if p(isig)<0.001    ptext = [ptext '*'];    end
            text((bar(ibar)+dif(isig))*1.1, (mu(isig,ibar))*1.1, ptext,'fontsize',20);
        end
    end
end




%interaction by subject
figure('color','white','name',blkstr);
hold on
myconds = {'dInd','dDep','dDepG'};
lconds = {'openLoop','closedLoop'};
mymin = 1; mymax = 0; a = [];
mycolors = mycolors(2:end,:);
for ilcond = 1:2
    a(ilcond) = subplot(1,2,ilcond);
    tmp = squeeze(avgdepdif(:,:,ilcond));
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
figure('color','white','name',blkstr);
mymin = min([min(min(nanmean(oProb,3))), min(min(nanmean(cProb,3)))]);
mymax = max([max(max(nanmean(oProb,3))), max(max(nanmean(cProb,3)))]);
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
imsc(nanmean(cProb,3), [mymin mymax], 'jet', [.5 .5 .5]); colorbar;
set(gca, 'DataAspectRatioMode', 'auto');
set(h,'ytick',1:1:size(cProb,2));
set(h,'yticklabel',probDO{1});
set(h,'xtick', 1:size(cProb,1));
set(h,'xticklabel',probDO{2});
title('ClosedLoop accuracy');
set(gca,'fontsize',16);




%plot by type
h = figure('color','white','name',blkstr);
mycolors = distinguishable_colors(8,{'w','k'});
hold on
mymin = 1;
a = [];
stops = [-1.5 -.5 .5 1.5];
for iaxis = 1:2
    a(iaxis) = subplot(1,2,iaxis);
    if iaxis == 1, idim = 2; else idim = 1; end %hack to fix mean dimenison mismatch
    omu = squeeze(nanmean(oProb,idim)); if size(omu,1)~= 4, omu = omu'; end
    cmu = squeeze(nanmean(cProb,idim)); if size(cmu,1)~=4, cmu = cmu'; end
    mu = cat(3,omu,cmu);
    mu = permute(mu, [3 1 2]);
    
    [xtick,hb,he]=errorbar_groups(mean(mu,3)', crit*ste(mu,3)','bar_names',{'OpenLoop','ClosedLoop'},'FigID', h, 'AxID', a(iaxis),...
        'bar_colors',mycolors(1+((iaxis-1)*size(mycolors,1)/2):(iaxis*size(mycolors,1)/2),:),...
        'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',5});
    
    legend(probDO{iaxis}');
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



if isfield(in,'ntacc') %is there non-targ accuracy? 
    nEvts = size(in.ntCueTarg,1);
    ol = in.ntCueTarg(1:nEvts/2,:,:);
    cl = in.ntCueTarg(nEvts/2+1:end,:,:);
    olacc = squeeze(mean(ol,1))';
    clacc = squeeze(mean(cl,1))';
    tmp = cat(3,olacc,clacc); tmp = permute(tmp,[3,1,2]);
    mycolors = distinguishable_colors(2,{'w','k'});
    [xtick,hb,he] = errorbar_groups(mean(tmp,3), crit*ste(tmp,3),'bar_names',in.ntCueTargDO{3},...
        'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',5});
    legend({'openLoop','closedLoop'});
    hold on
    stops = [-.5 .5];
    for igrp = 1:length(xtick)
        plot(repmat(xtick(igrp)+stops(1),size(tmp,3),2),squeeze(tmp(1,igrp,:)),'k.','markersize',20)
        plot(repmat(xtick(igrp)+stops(2),size(tmp,3),2),squeeze(tmp(2,igrp,:)),'k.','markersize',20)
    end

    ntacc = mean(in.ntacc,2)';
    ntste = ste(in.ntacc,2);
    for icond = 1:2
        shadedErrorBar(xlim, repmat(ntacc(icond),2,1), repmat(crit.*ntste(icond),2,1),...
            {'--','linewidth',2,'color',mycolors(icond,:),'markerfacecolor',mycolors(1,:)},1);
    end
    
    set(gca,'fontsize',20);
    ylabel('non-target accuracy');
    
end

if cfg.blkplts
    %testing block plots
    blkdiff = cat(3,squeeze(in.avgdepdif(:,1,1,:)),squeeze(in.avgdepdif(:,1,2,:)));
    if cfg.ci
        crit = tinv(.975,size(blkdiff,1)-1);
    else
        crit = 1;
    end
    h = figure('color','white');
    mycolors = get(gca,'defaultAxesColorOrder');
    errorbar_groups(squeeze(mean(blkdiff,1))',crit.*squeeze(ste(blkdiff,1))',...
        'bar_names',{'Blk1','Blk2','Blk3','Blk4'},'bar_colors',mycolors(5:end,:),'FigID', h);
    hold on
    if cfg.ci
        avgcrit = tinv(.975,4*size(blkdiff,1)-1);
    else
        avgcrit = 1;
    end
    for icond = 1:size(blkdiff,3)
        condmean =  mean(reshape(blkdiff(:,:,icond),numel(blkdiff(:,:,icond)),1));
        condste = ste(reshape(blkdiff(:,:,icond),numel(blkdiff(:,:,icond)),1));
        shadedErrorBar(xlim,repmat(condmean,1,2),repmat(avgcrit*condste,1,2),...
            {'--','linewidth',2,'color',mycolors(4+icond,:),'markerfacecolor',mycolors(1,:)},1);
    end
    ylabel('data-independent model');
    legend('openLoop','closedLoop');
    set(gca,'fontsize',20);
    
    olacc = deevResAcc(in.res(1:18,:,:,:))';
    clacc = deevResAcc(in.res(19:end,:,:,:))';
    oldep = squeeze(in.avgdep(1,1,:,:));
    cldep = squeeze(in.avgdep(2,1,:,:));
    
    mycolors = distinguishable_colors(4);
    blkstr = {'blk1','blk2','blk3','blk4'};
    h = figure('color','white');
    
    ax = subplot(2,2,1);
    errorbar_groups(mean(olacc),crit.*ste(olacc),'bar_colors',mycolors(1,:),'bar_names', blkstr, 'FigID', h, 'AxID', ax,...
        'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',5});
    ylabel('accuracy'); title('openloop'); ylim([0,1]); set(gca,'fontsize',20);
    
    ax = subplot(2,2,2);
    errorbar_groups(mean(clacc),crit.*ste(clacc),'bar_colors',mycolors(2,:),'bar_names', blkstr, 'FigID', h, 'AxID', ax,...
        'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',5});
    title('closedloop'); ylim([0,1]); set(gca,'fontsize',20);
    
    ax = subplot(2,2,3);
    errorbar_groups(mean(oldep),crit.*ste(oldep),'bar_colors',mycolors(3,:),'bar_names', blkstr, 'FigID', h, 'AxID', ax,...
        'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',5});
    ylabel('dependency'); ylim([.5,1]); set(gca,'fontsize',20);
    
    ax = subplot(2,2,4);
    errorbar_groups(mean(cldep),crit.*ste(cldep),'bar_colors',mycolors(4,:),'bar_names', blkstr, 'FigID', h, 'AxID', ax,...
        'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',5});
    ylim([.5,1]); set(gca,'fontsize',20);
   
end