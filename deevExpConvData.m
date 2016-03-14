function deevExpConvData(substr,cfg)
% exports data to be used in convolution net
%
% input:
%   substr: string, subject data to export
%   cfg: config struct
%       .adfile: string, analysis details file, def:'/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/v1/ft_data/OL_CL_encLoc_encPer_encObj_encAni_encNotLoc_encNotPer_encNotObj_encNotAni_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
%       .chansel: string, which channel group to use in export, def:'noEyeABH'
%       
% output: saves mat file with exported training and test data and
% corresponding labels in dirs.dataroot/dirs.saveDirStem/ft_exp/substr
%   

if ~exist('cfg','var'),          cfg = [];                              end
if ~isfield(cfg,'adfile'),       cfg.adFile = '/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/v1/ft_data/OL_CL_encLoc_encPer_encObj_encAni_encNotLoc_encNotPer_encNotObj_encNotAni_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat'; end
if ~isfield(cfg,'chansel'),      cfg.chansel = 'noEyeABH';              end

adFile = cfg.adFile; chansel = cfg.chansel;

%load analysis details
load(adFile) %contains exper, ana, dir structures

%select specific subject
badInd = zeros(size(exper.subjects));
badInd(~strcmp(exper.subjects,substr)) = 1;
exper = nk_rmSubs(exper,badInd);

%load data parameters
%cfg = [];

cfg.loadMethod = 'seg';
cfg.latency = 'all';
cfg.frequency = 'all';

cfg.keeptrials = 'yes';
cfg.equatetrials = 'no';

cfg.rmPreviousCfg = true;

% type of input (used in the filename to load)
cfg.ftype = 'pow';

% type of output: 'pow', 'coh', 'phase'
cfg.output = 'pow';

% transformation: 'log10', 'log', 'vec'
cfg.transform = '';

% normalization of single or average trials
% cfg.norm_trials = 'single'; % Grandchamp & Delorme (2011)
cfg.norm_trials = 'single';

% baseline type
% % 'zscore', 'absolute', 'relchange', 'relative', 'db'
 cfg.baseline_type = 'zscore';
% cfg.baseline_type = 'db';
% cfg.baseline_type = 'absolute';
% cfg.baseline_type = 'relchange';
% cfg.baseline_type = 'relative';

%cfg.baseline_time = [-0.3 -0.1];
cfg.baseline_time = [];

% at what data stage should it be baseline corrected?
cfg.baseline_data = 'pow';
% mod is not an option
% % cfg.baseline_data = 'mod';

%cfg.saveFile = true;
cfg.saveFile = false;

% only keep induced data by removing evoked?
cfg.rmevoked = 'no';
cfg.rmevokedfourier = 'no';
cfg.rmevokedpow = 'no';
% % baseline using all events
% cfg.baseline_events = 'all';

ana = mm_ft_elecGroups(ana);
chanstr = ana.elecGroups{strcmp(ana.elecGroupsStr,chansel)};

%%
%load encoding data
evtStr = 'encLoc';
ana.eventValues = {{evtStr}};
pos =  mm_ft_loadData_multiSes(cfg,exper,dirs,ana);
pos = pos.ses1.(evtStr).sub.data;

evtStr = 'encNotLoc';
ana.eventValues = {{evtStr}};
neg =  mm_ft_loadData_multiSes(cfg,exper,dirs,ana);
neg = neg.ses1.(evtStr).sub.data;

tsel = pos.time > -.1 & pos.time < pos.time(end);
t = pos.time(tsel); t = [t(1) t(end)];
fsel = pos.freq >= 4 & pos.freq <= 50;
f = pos.freq(fsel); f = [f(1) f(end)];

%prep for export
train = cat(1, pos.powspctrm(:,ismember(pos.label,chanstr),fsel,tsel),...
               neg.powspctrm(:,ismember(neg.label,chanstr),fsel,tsel));
if sum(isnan(train(:))) > 0 error('nan values in train data'); end
%reshape chan x freq x time x trials
%train = permute(train,[2,3,4,1]);
           
lbl = cat(1, pos.trialinfo, neg.trialinfo);
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
test_fsel = pos.freq >= 4 & pos.freq <= 50;

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

%% save out for python

savedir = fullfile(dirs.dataroot, dirs.saveDirStem,'ft_exp',substr);
if ~exist(savedir,'dir'), mkdir(savedir); end
fname = sprintf('%s_%s_freq%.02fto%.02f_time%.02fto%.02f.mat',substr,chansel,f(1),f(2),t(1),t(2));
fprintf('saving %s train data and labels...',substr);
save(fullfile(savedir,fname),'train','train_lbl','test','test_lbl','test_trialinfo','train_trialinfo','cfg');
fprintf('done\n\n');





