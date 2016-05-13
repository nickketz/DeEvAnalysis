function deevExpConvDataLong(substr,cfg)
% exports data to be used in convolution net
%
% input:
%   substr: string, subject data to export
%   cfg: config struct
%       .adfile: string, analysis details file, def:'/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/long/ft_data/rFix_OL_CL_rlOL_rlCL_Enc_eq0_art_ftAuto/tfr_wavelet_w4_fourier_3_50/analysisDetails.mat'
%       .chansel: string, which channel group to use in export, def:'noEyeABH'
%       .doRL: bool, export response locked testing data, otherwise do stim locked, def:0
%       .twin: array, start and stop (non-inclusive) time points taken from encoding data
%
% output: saves mat file with exported training and test data and
% corresponding labels in dirs.dataroot/dirs.saveDirStem/ft_exp/substr
%

if ~exist('cfg','var'),          cfg = [];                              end
if ~isfield(cfg,'adfile'),       cfg.adFile = '/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/long/ft_data/rFix_OL_CL_rlOL_rlCL_Enc_eq0_art_ftAuto/tfr_wavelet_w4_fourier_3_50/analysisDetails.mat'; end
if ~isfield(cfg,'chansel'),      cfg.chansel = 'noEyeABH';              end
if ~isfield(cfg,'doRL'),         cfg.doRL = 1;                          end
if ~isfield(cfg,'twin'),         cfg.twin = [.15,3.33];                  end

adFile = cfg.adFile; chansel = cfg.chansel; doRL = cfg.doRL; twin = cfg.twin;

%load analysis details
load(adFile) %contains exper, ana, dir structures

%select specific subject
badInd = zeros(size(exper.subjects));
badInd(~strcmp(exper.subjects,substr)) = 1;
exper = nk_rmSubs(exper,badInd);

%load data parameters
%cfg = [];

cfg = deevLoadConfig;
cfg.ftype = 'fourier';
cfg.keeptrials = 'yes';
cfg.baseline_type = 'zscore';
cfg.norm_trials = 'average';
cfg.baseline_time = [-.3 -.1];


ana = mm_ft_elecGroups(ana);
chanstr = ana.elecGroups{strcmp(ana.elecGroupsStr,chansel)};

%%
%load encoding data
evtStr = 'Enc';
ana.eventValues = {{evtStr}};
pos =  mm_ft_loadData_multiSes(cfg,exper,dirs,ana);
pos = pos.ses1.(evtStr).sub.data;

tsel = pos.time > twin(1) & pos.time < twin(2);
t = pos.time(tsel); t = [t(1) t(end)];
fsel = pos.freq >= 3 & pos.freq <= 48;
f = pos.freq(fsel); f = [f(1) f(end)];

%prep for export
train = pos.powspctrm(:,ismember(pos.label,chanstr),fsel,tsel);

if sum(isnan(train(:))) > 0 error('nan values in train data'); end
%reshape chan x freq x time x trials
%train = permute(train,[2,3,4,1]);

lbl = pos.trialinfo;
train_trialinfo = lbl;
lbl = lbl(:,8);%sttp column
lbl(:,2) = mod(lbl(:,1),10);
lbl(:,3) = (lbl(:,1) - lbl(:,2))/10;
train_lbl = zeros(size(lbl,1),4);
for i = 1:size(lbl,1)
    train_lbl(i,lbl(i,2)) = 1;
    train_lbl(i,lbl(i,3)) = 1;
end
clear pos neg

%%
if ~doRL
    
    %load testing data
    evtStr = 'CL';
    ana.eventValues = {{evtStr}};
    pos =  mm_ft_loadData_multiSes(cfg,exper,dirs,ana);
    pos = pos.ses1.(evtStr).sub.data;
    
    evtStr = 'OL';
    ana.eventValues = {{evtStr}};
    neg =  mm_ft_loadData_multiSes(cfg,exper,dirs,ana);
    neg = neg.ses1.(evtStr).sub.data;
    
    test_tsel = pos.time > -.1 & pos.time < pos.time(end);
    test_fsel = pos.freq >= f(1) & pos.freq <= f(2);
    
    if sum(test_fsel~=fsel)>0, error('mismatched frequency selection between training and test'); end
    
    %prep for export
    test = cat(1, pos.powspctrm(:,ismember(pos.label,chanstr),test_fsel,test_tsel),...
        neg.powspctrm(:,ismember(neg.label,chanstr),test_fsel,test_tsel));
    if sum(isnan(test(:)))>0, error('nan values in test data'); end
    %reshape as chans x freq x time x trials
    %test = permute(test,[2,3,4,1]);
    
    lbl = cat(1, pos.trialinfo, neg.trialinfo);
    test_trialinfo = lbl;
    lbl = lbl(:,8);%sttp column
    lbl(:,2) = mod(lbl(:,1),10);
    lbl(:,3) = (lbl(:,1) - lbl(:,2))/10;
    test_lbl = zeros(size(lbl,1),4);
    for i = 1:size(lbl,1)
        test_lbl(i,lbl(i,2)) = 1;
        test_lbl(i,lbl(i,3)) = 1;
    end
    clear pos neg
  
    
    %% load response locked test data
    
else %doRL
    
    ana = deevFixTrlOrd(ana);
    
    blana = [];
    blana.eventValues = {{'rlOL','OL'}};
    blana.eventValuesSplit = {{ {'rlOLrt'} {'rFix'} }};
    blana.trl_expr = {{ {'sttp==1'} ...
        {'sttp==1 | sttp==2'};
        }};
    blana.trl_order = ana.trl_order;
    inana(1) = blana;
    
    blana = [];
    blana.eventValues = {{'rlCL','CL'}};
    blana.eventValuesSplit = {{ {'rlCLrt'} {'rFix'} }};
    blana.trl_expr = {{
        {'sttp==2'} ...
        {'sttp==1 | sttp==2'};
        }};
    blana.trl_order = ana.trl_order;
    inana(2) = blana;
    
    cfg_bl = [];
    
    trgCond = {'rlOLrt','rlCLrt'};
    blCond = repmat({'rFix'},size(trgCond));
    
    cfg_bl.loadcfg = deevLoadConfig;
    cfg_bl.loadcfg.ftype = 'fourier';
    cfg_bl.loadcfg.loadMethod = 'trialinfo';
    cfg_bl.loadcfg.keeptrials = 'yes';
    cfg_bl.loadcfg.baseline_type = '';
    cfg_bl.loadcfg.norm_trial = 'average';
    cfg_bl.loadcfg.baseline_time = [];
    
    cfg_bl.keeptrials = 'yes';
    cfg_bl.method = 'zscore';
    cfg_bl.twin = [-.3 -.1];
    cfg_bl.param = 'powspctrm';
    
    cfg_bl.nullRT = true;
    cfg_bl.bufRT = .3;
    
    for icond = 1:length(trgCond)
        cfg_bl.colinds = true(size(ana.trl_order.rlCL));%they're all the same size anyway
        cfg_bl.colinds(1) = 0; %condn col can't be used in matching
        cfg_bl.blCond = blCond{icond};
        cfg_bl.trgCond = trgCond{icond};
        
        [bldata, exper] = deevTrialBL(exper,dirs,inana(icond),cfg_bl);
        tdata.([cfg_bl.trgCond]) = bldata;
    end
    pos = tdata.('rlCLrt').sub.data;
    neg = tdata.('rlOLrt').sub.data;
    clear tdata;
    
    test_tsel = pos.time >= -t(2) & pos.time <= -t(1);
    test_fsel = pos.freq >= f(1) & pos.freq <= f(2);
    test_time = pos.time;
    if sum(test_fsel~=fsel)>0, error('mismatched frequency selection between training and test'); end
    
    %prep for export
    test = cat(1, pos.powspctrm(:,ismember(pos.label,chanstr),test_fsel,test_tsel),...
        neg.powspctrm(:,ismember(neg.label,chanstr),test_fsel,test_tsel));
    %if sum(isnan(test(:)))>0, error('nan values in test data'); end
    %reshape as chans x freq x time x trials
    %test = permute(test,[2,3,4,1]);
    
    lbl = cat(1, pos.trialinfo, neg.trialinfo);
    test_trialinfo = lbl;
    icol = ismember(ana.trl_order.CL,'ordn');
    lbl = lbl(:,icol);%ordn column
    lbl(:,2) = mod(lbl(:,1),10);
    lbl(:,3) = (lbl(:,1) - lbl(:,2))/10;
    if any(lbl(:)<1), error('element decomposition is broken'); end
    test_lbl = zeros(size(lbl,1),4);
    for i = 1:size(lbl,1)
        test_lbl(i,lbl(i,2)) = 1;
        test_lbl(i,lbl(i,3)) = 1;
    end
    clear pos neg
    
    cfg.cfg_bl = cfg_bl;
    
end
if size(test,2)~=size(train,2)
    warning('mismatch between train and test sizes')
end
%% save out for python

savedir = fullfile(dirs.dataroot, dirs.saveDirStem,'ft_exp',substr);
if ~exist(savedir,'dir'), mkdir(savedir); end
fname = sprintf('%s_%s_freq%.02fto%.02f_time%.02fto%.02f.mat',substr,chansel,f(1),f(2),t(1),t(2));
fprintf('saving %s train data and labels...',substr);
save(fullfile(savedir,fname),'train','train_lbl','test','test_lbl','test_trialinfo','train_trialinfo','cfg');
fprintf('done\n\n');





