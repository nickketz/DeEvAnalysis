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
%
%
%

%set defaults
if ~exist('cfg','var')      cfg = [];                     end
if ~isfield(cfg,'ci')       cfg.ci = 1;                   end
if ~isfield(cfg,'blk')      cfg.blk = 'avg';              end


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
    
end
    
nsubs = size(in.logdata.data.Subj,2);

if cfg.ci
    crit = tinv(.975,nsubs-1);
else
    crit = 1;
end

barstr = {'data','indp','dpnd','dpnd+g'};

%plot dependency
h = figure('color','white','name',blkstr);
mycolors = get(gca,'defaultAxesColorOrder');
errorbar_groups(mean(avgdep,3)',crit*ste(avgdep,3)','bar_names',{'OpenLoop','ClosedLoop'},'bar_colors',mycolors,'FigID',h,...
    'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',5});
legend(barstr,'location','best');
title('Dependency frm log');
ylim([min(min(mean(avgdep,3)))-.1 1]);
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
h = figure('color','white','name',blkstr);
mycolors = get(gca,'defaultAxesColorOrder');
errorbar_groups(squeeze(mean(avgdepdif,1)),crit*squeeze(ste(avgdepdif,1)),'bar_names',{'OpenLoop','ClosedLoop'},'bar_colors',mycolors(2:end,:),'FigID',h,...
    'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',5});
legend({'dIndp','dDpnd','dDpnd+G'},'location','best');
title('Delta Dependency');
%ylim([min(depdif(:)) max(depdif(:))]);
set(gca,'fontsize',20);


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
