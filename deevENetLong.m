function [b,fit,p,l,auc] = deevENetLong(substr,cfg)
% exports data to be used in convolution net
%
% input:
%   substr: string, subject data to export
%   cfg: config struct
%       .adfile: string, analysis details file, def:'/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/long/ft_data/rFix_OL_CL_rlOL_rlCL_Enc_eq0_art_ftAuto/tfr_wavelet_w4_fourier_3_50/analysisDetails.mat'
%       .chansel: string, which channel group to use in export, def:'noEyeABH'
%
% output: saves mat file with exported training and test data and
% corresponding labels in dirs.dataroot/dirs.saveDirStem/ft_exp/substr
%

if ~exist('cfg','var'),          cfg = [];                              end
if ~isfield(cfg,'adfile'),       cfg.adFile = '/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/long/ft_data/rFix_OL_CL_rlOL_rlCL_Enc_eq0_art_ftAuto/tfr_wavelet_w4_fourier_3_50/analysisDetails.mat'; end
if ~isfield(cfg,'chansel'),      cfg.chansel = 'noEyeABH';              end
if ~isfield(cfg,'tlen'),         cfg.tlen = 5;                          end
if ~isfield(cfg,'tslide'),       cfg.tslide = cfg.tlen-1;               end
if ~isfield(cfg,'balance'),      cfg.balance = 1;                       end


adFile = cfg.adFile; chansel = cfg.chansel; tlen = cfg.tlen; tslide = cfg.tslide; balance = cfg.balance;

%load analysis details
load(adFile) %contains exper, ana, dir structures

%select specific subject
badInd = zeros(size(exper.subjects));
badInd(~strcmp(exper.subjects,substr)) = 1;
exper = nk_rmSubs(exper,badInd);

%load data parameters
%cfg = [];

cfg_load = deevLoadConfig;
cfg_load.ftype = 'fourier';
cfg_load.keeptrials = 'yes';
cfg_load.baseline_type = '';
cfg_load.norm_trial = 'average';
cfg_load.baseline_time = [-.3 -.1];

ana = mm_ft_elecGroups(ana);
ana = deevFixTrlOrd(ana);
chanstr = ana.elecGroups{strcmp(ana.elecGroupsStr,chansel)};


%%
%load encoding data
evtStr = 'Enc';
ana.eventValues = {{evtStr}};
pos =  mm_ft_loadData_multiSes(cfg_load,exper,dirs,ana);
pos = pos.ses1.(evtStr).sub.data;

tsel = pos.time > -.1 & pos.time < pos.time(end);
t = pos.time(tsel); t = [t(1) t(end)];
fsel = pos.freq >= 4 & pos.freq <= 50;
f = pos.freq(fsel); f = [f(1) f(end)];

%prep for export
train = pos.powspctrm(:,ismember(pos.label,chanstr),fsel,tsel);

if sum(isnan(train(:))) > 0 error('nan values in train data'); end
%reshape chan x freq x time x trials
%train = permute(train,[2,3,4,1]);

lbl = pos.trialinfo;
train_trialinfo = lbl;
icol = ismember(ana.trl_order.CL,'ordn');
lbl = lbl(:,icol);%ordn column
lbl(:,2) = mod(lbl(:,1),10);
lbl(:,3) = (lbl(:,1) - lbl(:,2))/10;
if any(lbl(:)<1) error('element decomposition is broken'); end
train_lbl = zeros(size(lbl,1),4);
for i = 1:size(lbl,1)
    train_lbl(i,lbl(i,2)) = 1;
    train_lbl(i,lbl(i,3)) = 1;
end
time = pos.time;
clear pos neg

%%
% %load testing data
% evtStr = 'CL';
% ana.eventValues = {{evtStr}};
% pos =  mm_ft_loadData_multiSes(cfg,exper,dirs,ana);
% pos = pos.ses1.(evtStr).sub.data;
% 
% evtStr = 'OL';
% ana.eventValues = {{evtStr}};
% neg =  mm_ft_loadData_multiSes(cfg,exper,dirs,ana);
% neg = neg.ses1.(evtStr).sub.data;
% 
% test_tsel = pos.time > -.1 & pos.time < pos.time(end);
% test_fsel = pos.freq >= 4 & pos.freq <= 50;
% 
% if sum(test_fsel~=fsel)>0, error('mismatched frequency selection between training and test'); end
% 
% %prep for export
% test = cat(1, pos.powspctrm(:,ismember(pos.label,chanstr),test_fsel,test_tsel),...
%     neg.powspctrm(:,ismember(neg.label,chanstr),test_fsel,test_tsel));
% if sum(isnan(test(:)))>0, error('nan values in test data'); end
% %reshape as chans x freq x time x trials
% %test = permute(test,[2,3,4,1]);
% 
% lbl = cat(1, pos.trialinfo, neg.trialinfo);
% test_trialinfo = lbl;
% lbl = lbl(:,8);%ordn column
% lbl(:,2) = mod(lbl(:,1),10);
% lbl(:,3) = (lbl(:,1) - lbl(:,2))/10;
% if any(lbl(:)<1), error('element decomposition is broken'); end
% test_lbl = zeros(size(lbl,1),4);
% for i = 1:size(lbl,1)
%     test_lbl(i,lbl(i,2)) = 1;
%     test_lbl(i,lbl(i,3)) = 1;
% end
% clear pos neg

%% load response locked testing data
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

test_tsel = pos.time > -3 & pos.time < .1;
test_fsel = pos.freq >= 4 & pos.freq <= 50;
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
 

%% run elasticnet
%parpool();
tstart = 1:tslide:size(train,4)-tlen;
tstart_test = 1:1:size(test,4)-tlen;

conds = {'encLoc','encPer','encObj','encAni'};
for icond = 1:length(conds)

    %balance training data?
    if balance
        if ~exist('otrain','var')
            otrain = train;
            otrain_lbl = train_lbl;
        else
            train = otrain;
            train_lbl = otrain_lbl;
        end
        tmp = train_lbl(:,icond);
        if mean(tmp)>.5 %resample neg 
            resampclass = 0;
        else %resample pos
            resampclass = 1;
        end
        ind = find(tmp==resampclass);
        nsamp = sum(tmp~=resampclass) - sum(tmp==resampclass);
        samps = ind(randi(length(ind),nsamp,1));
        fprintf('\nbalancing conditions by resampling %d trials from class %d',nsamp,resampclass);
        train = cat(1,train,train(samps,:,:,:));
        train_lbl = cat(1,train_lbl,train_lbl(samps,:));
        if mean(train_lbl(:,icond)) ~= .5, error('balancing training data failed'); end
    end
    
    fprintf('\n\n\ncondition: %s',conds{icond});
    b = cell(length(tstart),1);
    fit = cell(length(tstart),1);
    y = train_lbl(:,icond);
    
    l = test_lbl(:,icond);
    maxp = nan(length(l),length(tstart));
    maxt = nan(length(l),length(tstart));
    meanp = maxp;
    maxauc = nan(length(tstart),1);
    meanauc = maxauc;
    maxdp = maxauc;
    meandp = maxauc;
    allp = [];
                
    for i = 1:length(tstart)
        tind = tstart(i):tstart(i)+tlen;
        fprintf('\ntraining on time window: %.02f to %.02f (%i/%i)',time(tind(1)),time(tind(end)),i,length(tstart));
        x = reshape(train(:,:,:,tind),[size(train,1),numel(train(1,:,:,tind))]);
        x = zscore(x);
        opts = statset('UseParallel',true);
        [b{i},fit{i}] = lassoglm(x,y,'binomial','alpha',.9,'Options',opts,'NumLambda',10,'CV',10);
        [~,icol] = min(fit{i}.Deviance);
        trainp = 1./(1+exp(-(fit{i}.Intercept(icol)+b{i}(:,icol)'*x')));
        tp = trainp*y; fp = trainp*(1-y); fn = (1-trainp)*y; tn = (1-trainp)*(1-y);
        trainerr = 1-(tp+tn)/(tp+tn+fp+fn);
        fprintf('\ntrain error: binary:%.2f norm:%.04f',mean(abs(round(trainp)'-y)),trainerr);

        p = nan(length(l),length(tstart_test));
        for j = 1:length(tstart_test)
            jtind = tstart_test(j):tstart_test(j)+tlen;
            x = reshape(test(:,:,:,jtind),[size(test,1),numel(test(1,:,:,jtind))]);
            x = nanzscore(x);
            p(:,j) = 1./(1+exp(-(fit{i}.Intercept(icol)+b{i}(:,icol)'*x')));
        end
        allp = cat(3,allp, p);
        
        %do max diff       
        [~,maxt(:,i)] = nanmax(abs(p-(1-p)),[],2);
        for ip = 1:size(maxt,1), maxp(ip,i) = p(ip,maxt(ip,i)); end
        [~,~,~,maxauc(i)] = perfcurve(l,maxp(:,i),1);
        binvar = maxp(:,i);%round(maxp(:,i));
        hit = binvar'*l;%sum(binvar==1 & l==1)/sum(l==1);
        fa = (1-binvar)'*l;%sum(binvar==1 & l==0)/sum(l==0);
        fn = binvar'*(1-l);
        tn = (1-binvar)'*(1-l);
        err = 1-(hit+tn)/(hit+tn+fp+fn);%mean(abs(l-binvar));
        maxdp(i) = norminv(sum((binvar.*l)>.5)/sum(l)) - norminv(sum((binvar.*(1-l))>.5)/sum(1-l));
        fprintf('\nmaxeval: auc %.02f; err %.02f; prec %.02f; recall %.02f; dprime %.02f ',...
            maxauc(i),err,hit/(hit+fa),hit/(hit+fn),maxdp(i));
        maxt(:,i) = test_time(maxt(:,i));
        
        %do mean
        meanp(:,i) = nanmean(p,2);
        [~,~,~,meanauc(i)] = perfcurve(l,meanp(:,i),1);
        binvar = meanp(:,i);%round(meanp(:,i));
        hit = binvar'*l;%sum(binvar==1 & l==1)/sum(l==1);
        fa = (1-binvar)'*l;%sum(binvar==1 & l==0)/sum(l==0);
        fn = binvar'*(1-l);
        tn = (1-binvar)'*(1-l);
        err = 1-(hit+tn)/(hit+tn+fp+fn);%mean(abs(l-binvar));
        meandp(i) = norminv(sum((binvar.*l)>.5)/sum(l)) - norminv(sum((binvar.*(1-l))>.5)/sum(1-l));
        fprintf('\nmeaneval: auc %.02f; err %.02f; hit %.02f; fa %.02f; dprime %.02f ',...
            meanauc(i),err,hit/(hit+fa),hit/(hit+fn),meandp(i));
        
        
    end

    %save    
    savedir = fullfile(dirs.saveDirProc,substr,ana.sessionNames{1});
    if ~exist(savedir,'dir'), mkdir(savedir); end
    fname = sprintf('Enet_%s_%dtlen_%dtslide.mat',conds{icond},tlen,tslide);
    fprintf('\n\nsaving to file:\n %s',fname);
    save(fullfile(savedir,fname),'train_lbl','test_lbl','test_trialinfo','train_trialinfo','cfg',...
        'f','t','chansel','maxp','maxt','meanp','allp','meanauc','maxauc','meandp','maxdp','fit','b');
end



