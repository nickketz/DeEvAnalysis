function deevRSA(substr,cfg)
% preliminary representation similarity analysis (RSA)
%
% input:
%   substr: string, subject data to export
%   cfg: config struct
%       .adfile: string, analysis details file, def:'/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/long/ft_data/rFix_OL_CL_rlOL_rlCL_Enc_eq0_art_ftAuto/tfr_wavelet_w4_fourier_3_50/analysisDetails.mat'
%       .chansel: string, which channel group to use in export, def:'noEyeABH'
%       .time: array, time window to do analysis on, def: [-.5, -.3]
%       .freq: array, freq window to do analysis on, def: [4,50]
%       .doRL: bool, export response locked testing data, otherwise do stim locked, def:1
%
% output: saves mat file with RSA matrix
%

if ~exist('cfg','var'),          cfg = [];                              end
if ~isfield(cfg,'adfile'),       cfg.adFile = '/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/long/ft_data/rFix_OL_CL_rlOL_rlCL_Enc_eq0_art_ftAuto/tfr_wavelet_w4_fourier_3_50/analysisDetails.mat'; end
if ~isfield(cfg,'chansel'),      cfg.chansel = 'noEyeABH';              end
if ~isfield(cfg,'doRL'),         cfg.doRL = 0;                          end
if ~isfield(cfg,'time'),         cfg.time = [0.48 1.32];                 end
if ~isfield(cfg,'freq'),         cfg.freq = [4 8];                     end

adFile = cfg.adFile; chansel = cfg.chansel; doRL = cfg.doRL; 
twin = cfg.time; fwin = cfg.freq;
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
cfg_load.norm_trials = 'average';
cfg_load.baseline_time = [-.3 -.1];


ana = mm_ft_elecGroups(ana);
chanstr = ana.elecGroups{strcmp(ana.elecGroupsStr,chansel)};
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

    test_tsel = pos.time > twin(1) & pos.time < twin(2);
    test_fsel = pos.freq >= fwin(1) & pos.freq <= fwin(2);
    test_time = pos.time(test_tsel);
    test_freq = pos.freq(test_fsel);
    
    %prep for export
    test = cat(1, pos.powspctrm(:,ismember(pos.label,chanstr),test_fsel,test_tsel),...
        neg.powspctrm(:,ismember(neg.label,chanstr),test_fsel,test_tsel));
      
    
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
    
    test_tsel = pos.time > twin(1) & pos.time < twin(2);
    test_fsel = pos.freq >= fwin(1) & pos.freq <= fwin(2);
    test_time = pos.time(test_tsel);
    test_freq = pos.freq(test_fsel);
    
    %prep for export
    test = cat(1, pos.powspctrm(:,ismember(pos.label,chanstr),test_fsel,test_tsel),...
        neg.powspctrm(:,ismember(neg.label,chanstr),test_fsel,test_tsel));
    
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
    clear pos neg
    
    saveStr = 'RL';
        
end

%% do corr 

scol = strcmp(ana.trl_order.OL,'trln');
bcol = strcmp(ana.trl_order.OL,'exbk');
stcol = strcmp(ana.trl_order.OL,'sttp');

% 
% %across condition similarity
% stim = [12, 21];
% dCorr = cell(length(stim), length(stim));
% sCorr = dCorr;
% dTrlInf = dCorr;
% fprintf('running across condition RSA\n');
% for iblock = 1:2
%     for istim = 1:length(stim)
%         %sub select for loc-per trials
%         tind = lbl(:,1) == stim(istim) & test_trialinfo(:,bcol)==iblock;
%         data = test(tind,:,:,:);
%         ti = test_trialinfo(tind,:);
%         if sum(tind)>1
%             [~,sind] = sort(ti(:,scol));
%             data = reshape(data(sind,:,:,:),size(data,1),numel(data(1,:,:,:)));
%             ti = ti(sind,:);
%             dTrlInf{istim,iblock} = ti;
%             dcorr = corrcoef(data','rows','pairwise'); %remove nans
%             scorr = sum(dcorr,2)-1;
%             %        [h,p,ci,stat] = ttest2(scorr(ti(:,stcol)==1),scorr(ti(:,stcol)==2));
%             dCorr{istim,iblock} = dcorr;
%             sCorr{istim,iblock} = scorr;
%         else
%             dCorr{istim,iblock} = nan;
%             sCorr{istim,iblock} = nan;
%             dTrlInf{istim,iblock} = nan;
%         end
%         %figure
%         %imagesc(dcorr);
%         %title(['stim ' num2str(stim(istim)) ' - blk ' num2str(iblock)]);
%     end
% end


%within events similarity
evtCorr = nan(36,2);
evtInds = nan(size(test_trialinfo,1),36,2);
evtCnt = evtCorr;
fprintf('running within event RSA\n');
for ievt = unique(test_trialinfo(:,scol))'
    for iblk = unique(test_trialinfo(:,bcol))'
        %gather trials from ievt and iblk
        ind = test_trialinfo(:,scol) == ievt & test_trialinfo(:,bcol) == iblk;
        if sum(ind)>1
            data = reshape(test(ind,:,:,:),sum(ind),numel(test(1,:,:,:)));
            evtInds(:,ievt,iblk) = ind;
            evtCnt(ievt,iblk) = sum(ind);
            tmp = tril(corrcoef(data','rows','pairwise'),-1);           
            evtCorr(ievt,iblk) = mean(sum(atanh(tmp(tmp~=0)),2));
            %tmp = tril(corrcoef(data','rows','pairwise'),-1);
            %evtCorr(ievt,iblk) = mean(tmp(tmp~=0));
        end
    end
end

        

%% save out for python

savedir = fullfile(dirs.saveDirProc,substr,exper.sesStr{1});
if ~exist(savedir,'dir'), mkdir(savedir); end
fname = sprintf('RSA_%s_%s_freq%.02fto%.02f_time%.02fto%.02f.mat',saveStr,chansel,fwin(1),fwin(2),twin(1),twin(2));
fprintf('saving %s RSA results...',substr);
save(fullfile(savedir,fname),'test','test_lbl','test_time','test_freq','test_trialinfo','cfg_in','cfg_load',...
                             'evtCorr','evtInds','evtCnt');%,...
%                             'dCorr','sCorr','dTrlInf');
fprintf('done\n\n');





