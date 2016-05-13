function [pred,orig,auc_v,auc_t] =  deevLoadConvData(substr,cfg)
% loads data exported from convolution net
%
% input:
%   substr: string, subject data to export
%   cfg: config struct
%       .adfile: string, analysis details file, def:'/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/v1/ft_data/OL_CL_encLoc_encPer_encObj_encAni_encNotLoc_encNotPer_encNotObj_encNotAni_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
%       .CNNdir: path to CNN data def:'~/work/DeEv_EEG/analysis/CNN'
%
% output: saves mat file with exported training and test data and
% corresponding labels in dirs.dataroot/dirs.saveDirStem/ft_exp/substr
%

if ~exist('cfg','var'),          cfg = [];                              end
if ~isfield(cfg,'adfile'),       cfg.adFile = '/home/nike3851/work/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/v1/ft_data/OL_CL_encLoc_encPer_encObj_encAni_encNotLoc_encNotPer_encNotObj_encNotAni_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat'; end
if ~isfield(cfg,'chansel'),      cfg.chansel = 'noEyeABH';              end
if ~isfield(cfg,'CNNdir'),       cfg.CNNdir = '/home/nike3851/work/DeEv_EEG/analysis/CNN'; end
if ~isfield(cfg,'doplots'),      cfg.doplots = 0;                       end

adFile = cfg.adFile; chansel = cfg.chansel; CNNdir = cfg.CNNdir;

%load analysis details
load(adFile) %contains exper, ana, dir structures

%select specific subject
badInd = zeros(size(exper.subjects));
badInd(~strcmp(exper.subjects,substr)) = 1;
exper = nk_rmSubs(exper,badInd);

%% load exported data to verify output matches
%f = [4.28 49.94]; t=[-0.08 1.96];

savedir = fullfile(CNNdir,substr(6:end),'Backup','test_pred.mat');
fprintf('loading %s predictions...',substr);
pred = load(fullfile(savedir),'test_p','train_p','datafile');
fprintf('done\n\n');

fprintf('loading %s train data and labels...',substr);
orig = load(pred.datafile,'train_lbl','test_lbl','train_trialinfo','test_trialinfo');
fprintf('done\n\n');
%
% tp = pred.train_p;
% tp = abs(min(tp(:))) + tp;
% tp = tp./max(tp(:));
% pred.train_p = tp;
%
% tp = pred.test_p;
% tp = abs(min(tp(:))) + tp;
% tp = tp./max(tp(:));
% pred.test_p = tp;

if cfg.doplots
    figure('color','white')
    c = get(gca,'colororder');
    
    plot([0 1], [0 1], '--k');
    hold on
end

%     [x,y,~,auc] = perfcurve(orig.train_lbl(:),pred.train_p(:),1);
%     h1=plot(x,y,'-b','linewidth',2);
%     text(x(round(length(x)/4)),y(round(length(y)/4)),['auc=' num2str(auc,'%.02f')],'fontsize',20,'color','blue');
%     [x,y,~,auc] = perfcurve(orig.test_lbl(:),pred.test_p(:),1);
%     h2=plot(x,y,'-g','linewidth',2);
%     text(x(round(length(x)/2)),y(round(length(y)/2)),['auc=' num2str(auc,'%.02f')],'fontsize',20,'color','green');
%
conds = {'loc','per','obj','ani'};
lconds = {};%{'train','val'};
cndcnt = 0;
auc_t = nan(length(conds),1); auc_v = auc_t;
for icond = 1:size(pred.test_p,2)
    cndcnt = cndcnt +1;
    [x,y,~,auc_t(icond)] = perfcurve(orig.train_lbl(:,icond),pred.train_p(:,icond),1);
    if cfg.doplots
        hc(cndcnt) = plot(x,y,'--','linewidth',2,'color',c(icond,:),'linewidth',2);
    end
    %text(x(round(length(x)/4)),y(round(length(y)/4)),['auc=' num2str(auc,'%.02f')],'fontsize',10,'color',c(icond,:));
    lconds = cat(2,lconds,[conds{icond} '_t']);
    cndcnt = cndcnt+1;
    [x,y,~,auc_v(icond)] = perfcurve(orig.test_lbl(:,icond),pred.test_p(:,icond),1);
    if cfg.doplots
        hc(cndcnt) = plot(x,y,':','linewidth',2,'color',c(icond,:),'linewidth',2);
    end
    %text(x(round(length(x)/2)),y(round(length(y)/2)),['auc=' num2str(auc,'%.02f')],'fontsize',10,'color',c(icond,:));
    lconds = cat(2,lconds,[conds{icond} '_v']);
end
if cfg.doplots
    box off
    title(sprintf('subj %s performance',substr(end-1:end)));
    set(gca,'fontsize',20);
    %legend([h1,h2,hc],lconds,'fontsize',12);
    legend(hc,lconds,'fontsize',12);
end

end





