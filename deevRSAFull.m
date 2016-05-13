function deevRSAFull(substr,cfg)
% preliminary representation similarity analysis (RSA)
%
% input:
%   substr: string, subject data to export
%   cfg: config struct
%       .adfile: string, analysis details file, def:'/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/long/ft_data/rFix_OL_CL_rlOL_rlCL_Enc_eq0_art_ftAuto/tfr_wavelet_w4_fourier_3_50/analysisDetails.mat'
%       .doRL: bool, export response locked testing data, otherwise do stim locked, def:1
%
% output: saves mat file with RSA matrix
%

if ~exist('cfg','var'),          cfg = [];                              end
if ~isfield(cfg,'adfile'),       cfg.adFile = '/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/long/ft_data/rFix_OL_CL_rlOL_rlCL_Enc_eq0_art_ftAuto/tfr_wavelet_w4_fourier_3_50/analysisDetails.mat'; end
if ~isfield(cfg,'doRL'),         cfg.doRL = 0;                          end

adFile = cfg.adFile;
doRL = cfg.doRL; 
cfg_in = cfg;

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
cfg_load.baseline_type = 'zscore';
cfg_load.norm_trials = 'average';
cfg_load.baseline_time = [-.3 -.1];


ana = deevFixTrlOrd(ana);
%%
if ~doRL
    
    
    
    %load testing data
    evtStr = 'CL';
    rtcol = strcmp('rtmm',ana.trl_order.(evtStr));
    ana.eventValues = {{evtStr}};
    pos =  mm_ft_loadData_multiSes(cfg_load,exper,dirs,ana);
    pos = pos.ses1.(evtStr).sub.data;
    
    evtStr = 'OL';
    ana.eventValues = {{evtStr}};
    neg =  mm_ft_loadData_multiSes(cfg_load,exper,dirs,ana);
    neg = neg.ses1.(evtStr).sub.data;

    test_tsel = true(size(pos.time));
    test_fsel = true(size(pos.freq));
    test_time = pos.time(test_tsel);
    test_freq = pos.freq(test_fsel);
    
    %prep for export
    test = cat(1, pos.powspctrm(:,:,test_fsel,test_tsel),...
        neg.powspctrm(:,:,test_fsel,test_tsel));
      
    
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
    template = rmfield(pos,{'powspctrm','trialinfo','cumtapcnt'});
    clear pos neg
    
    %nan post response
    for itrl = 1:size(test_trialinfo)
        rt = test_trialinfo(itrl,rtcol);
        rtsel = test_time>=rt;
        if sum(rtsel)>0
            test(itrl,:,:,rtsel) = nan;
        end
    end
    
    saveStr = 'SL';
  
    
    %% load response locked test data
    
else %doRL
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
    
    cfg_load = [];
    
    trgCond = {'rlOLrt','rlCLrt'};
    blCond = repmat({'rFix'},size(trgCond));
    
    cfg_load.loadcfg = deevLoadConfig;
    cfg_load.loadcfg.ftype = 'fourier';
    cfg_load.loadcfg.loadMethod = 'trialinfo';
    cfg_load.loadcfg.keeptrials = 'yes';
    cfg_load.loadcfg.baseline_type = '';
    cfg_load.loadcfg.norm_trial = 'average';
    cfg_load.loadcfg.baseline_time = [];
    
    cfg_load.keeptrials = 'yes';
    cfg_load.method = 'zscore';
    cfg_load.twin = [-.3 -.1];
    cfg_load.param = 'powspctrm';
    
    cfg_load.nullRT = true;
    cfg_load.bufRT = .3;
    
    for icond = 1:length(trgCond)
        cfg_load.colinds = true(size(ana.trl_order.rlCL));%they're all the same size anyway
        cfg_load.colinds(1) = 0; %condn col can't be used in matching
        cfg_load.blCond = blCond{icond};
        cfg_load.trgCond = trgCond{icond};
        
        [bldata, exper] = deevTrialBL(exper,dirs,inana(icond),cfg_load);
        tdata.([cfg_load.trgCond]) = bldata;
    end
    pos = tdata.('rlCLrt').sub.data;
    neg = tdata.('rlOLrt').sub.data;
    clear tdata;
    
    test_tsel = true(size(pos.time));
    test_fsel = true(size(pos.freq));
    test_time = pos.time(test_tsel);
    test_freq = pos.freq(test_fsel);
    
    %prep for export
    test = cat(1, pos.powspctrm(:,:,test_fsel,test_tsel),...
        neg.powspctrm(:,:,test_fsel,test_tsel));
    
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
        test_lbl(i,lbl(i,3)) = 2;
    end
    template = rmfield(pos,{'powspctrm','trialinfo','cumtapcnt'});
    clear pos neg
    
    saveStr = 'RL';
        
end

%% do corr 

scol = strcmp(ana.trl_order.OL,'trln');
bcol = strcmp(ana.trl_order.OL,'exbk');
stcol = strcmp(ana.trl_order.OL,'sttp');


%within events similarity
evtCorr = nan(36,2);
evtInds = nan(size(test_trialinfo,1),36,2);
evtCnt = evtCorr;
fprintf('running within event RSA\n');
tsize = size(test);
data = nan([36,2,tsize(2:end)]);
for ievt = unique(test_trialinfo(:,scol))'
    for iblk = unique(test_trialinfo(:,bcol))'
        %gather trials from ievt and iblk
        ind = test_trialinfo(:,scol) == ievt & test_trialinfo(:,bcol) == iblk;
        if sum(ind)>2
            %variance across retrieval trials from a given event/blk
            data(ievt,iblk,:,:,:) = nanvar(test(ind,:,:,:),1); %need normalized variance?
        end
    end
end
%combine across blocks and OL/CL events
varOL = squeeze(reshape(data(1:18,:,:,:,:),[18*2,tsize(2:end)]));
cntOL = sum(sum(sum(sum(isnan(varOL),4),3),2) ~= numel(varOL(1,:)));
varCL = squeeze(reshape(data(19:end,:,:,:),[18*2,tsize(2:end)]));
cntCL = sum(sum(sum(sum(isnan(varCL),4),3),2) ~= numel(varCL(1,:)));


%% save out for python

savedir = fullfile(dirs.saveDirProc,substr,exper.sesStr{1});
if ~exist(savedir,'dir'), mkdir(savedir); end
fname = sprintf('RSAFull_%s.mat',saveStr);
fprintf('saving %s RSA results...',substr);
template.dimord = 'chan_freq_time';
save(fullfile(savedir,fname),'cfg_in','cfg_load',...
                             'varOL','varCL','cntOL','cntCL','template');
fprintf('done\n\n');

