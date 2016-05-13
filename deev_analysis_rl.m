%preliminary power analysis to check on flckr condtions
clear
%load analysis details
%short pow
%adFile = '/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/v1/ft_data/OL_CL_encLoc_encPer_encObj_encAni_encNotLoc_encNotPer_encNotObj_encNotAni_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
%long pow
%adFile = '/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/long/ft_data/OL_CL_encLoc_encPer_encObj_encAni_encNotLoc_encNotPer_encNotObj_encNotAni_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
%fourier long
adFile = '/home/nike3851/work/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/long/ft_data/rFix_OL_CL_rlOL_rlCL_Enc_eq0_art_ftAuto/tfr_wavelet_w4_fourier_3_50/analysisDetails.mat';

load(adFile)

%
%remove subs?
conds = {'rlOL','rlCL','OL','CL'};%ana.eventValues{1};
tmpexper = exper;
%find subs with 0 trials in any condition
nTrials = [];
mintrials = 20;
for icond = 1:length(conds)
    nTrials = cat(2,nTrials,exper.nTrials.(conds{icond}));
end
badInd = sum(nTrials<mintrials,2)>0;
%also remove sub20 for no behav data
%badInd(ismember(exper.subjects,{'DeEv_sub23','DeEv_sub24','DeEv_sub25'})) = 1;
if sum(badInd)>0
    bsubs = sprintf('%s, ',exper.subjects{badInd'});
    fprintf('\nRemoving subs: %s\n',bsubs(1:end-2));
   % exper = nk_rmSubs(exper,badInd);
else
    fprintf('\nall subs have atleast %d trials!\n',mintrials);
end


%make sure ana.trl_order is set correctly
ana = deevFixTrlOrd(ana);

exper.badSub = zeros(length(exper.subjects),1);
ana = mm_ft_elecGroups(ana);


%% load events - regular
ana.eventValues = { {'OL','CL'} };

cfg = deevLoadConfig;
cfg.keeptrials = 'no';
cfg.ftype = 'fourier';
cfg.norm_trials = 'average';
cfg.baseline_time = [-.3 -.1];

[data_pow,exper] =  mm_ft_loadData_multiSes(cfg,exper,dirs,ana);
exper.badSub = zeros(length(exper.subjects),1);
ana = mm_ft_elecGroups(ana);

     
%% events of interest - trialinfo
rtthres = 2.5;
ana.eventValues = {{'rFix','rFix',...
                    'rlOL','rlCL','OL','CL',...
                    'rlOL','rlCL','OL','CL',...
                    'rlOL','rlCL','OL','CL',...
                    'rlOL','rlCL','OL','CL'}};

ana.eventValuesSplit = {{ {'rFix_OL','rFix_CL'} {'rFix'}... 
                          {'rlOLrtcrt'} {'rlCLrtcrt'} {'OLrtcrt'} {'CLrtcrt'}...
                          {'rlOLcrt'} {'rlCLcrt'} {'OLcrt'} {'CLcrt'}...
                          {'rlOLrt'} {'rlCLrt'} {'OLrt'} {'CLrt'}...
                          {'rlOL'} {'rlCL'} {'OL'} {'CL'} }};
                      
ana.trl_expr = {{ {'sttp==1', 'sttp==2'} {'sttp==1 | sttp==2'} ... %{'rFix_OL','rFix_CL'} {'rFix'}
                  {['sttp==1 & crct==1 & rtmm>' num2str(rtthres)]} ... %{'rlOLrtcrt'}
                   {['sttp==2 & crct==1 & rtmm>' num2str(rtthres)]} ... %{'rlCLrtcrt'}  
                   {['sttp==1 & crct==1 & rtmm>' num2str(rtthres)]} ... %{'OLrtcrt'}
                   {['sttp==2 & crct==1 & rtmm>' num2str(rtthres)]} ... %{'CLrtcrt'}
                  {'sttp==1 & crct==1'} {'sttp==2 & crct==1'} ... %{'rlOLcrct'} {'rlCLcrct'}
                   {'sttp==1 & crct==1'} {'sttp==2 & crct==1'} ...%{'OLcrct'} {'CLcrct'}
                  {['sttp==1 & rtmm>' num2str(rtthres)]} ...%{'rlOLrt}
                   {['sttp==2 & rtmm>' num2str(rtthres)]}...%{'rlCLrt}
                   {['sttp==1 & rtmm>' num2str(rtthres)]} ...%{'OLrt}
                   {['sttp==2 & rtmm>' num2str(rtthres)]}...%{'CLrt}
                  {'sttp==1'} {'sttp==2'} {'sttp==1'} {'sttp==2'} }}; %{'rlOL'} {'rlCL'} {'OL'} {'CL'}




%% load the data
%define trials of interest
%ana.eventValues = {{'OL','CL'}};
%ana.eventValues = {exper.eventValues};

cfg = deevLoadzConfig;
cfg.keeptrials = 'no';
cfg.ftype = 'fourier';
cfg.baseline_time = [-.3 -.1];
cfg.loadMethod = 'trialinfo';

[data_pow,exper] =  mm_ft_loadData_multiSes(cfg,exper,dirs,ana);
exper.badSub = zeros(length(exper.subjects),1);
ana = mm_ft_elecGroups(ana);

tmp = {};
for isplit = 1:length(ana.eventValuesSplit{1})
    tmp = cat(2,tmp,ana.eventValuesSplit{1}{isplit});
end
ana.eventValues = {tmp};


%% condition based baseline
cfg_bl = [];
cfg_bl.twin = [0 .3];
trgCond = {'rlCL','rlOL','rlCLrt','rlOLrt','CL','OL',...
           'CLrt','OLrt','CLrtcrt','OLrtcrt','rlCLrtcrt','rlOLrtcrt'...
           'rlCLcrt','rlOLcrt','CLcrt','OLcrt'};
       
blCond = repmat({'rFix'},size(trgCond));

cfg_bl.method = 'zscore';
cfg_bl.param = 'powspctrm';
for icond = 1:length(trgCond)
    cfg_bl.blCond = blCond{icond};
    cfg_bl.trgCond = trgCond{icond};
    data_pow = deevCondBL(data_pow,cfg_bl);
end
%add new cond to eventValues
ana.eventValues = {fieldnames(data_pow.ses1)};
%select bl corrected data
ind = ~cellfun(@isempty,regexp(ana.eventValues{1},'^bl.*'));
ana.eventValues = {ana.eventValues{1}(ind)};

%% load trial based baseline correction
%make sure ana is setup correct
ana = deevFixTrlOrd(ana);
clear inana

blana = [];
blana.eventValues = {{'rlOL','OL'}};
blana.eventValuesSplit = {{ {'rlOLrtcrt'} {'rFix'} }};
blana.trl_expr = {{ {'sttp==1 & crct==1'} ... %{'rlOLrtcrt'}
                    {'sttp==1 | sttp==2'}; 
                 }};     
blana.trl_order = ana.trl_order;
inana{1} = blana;

blana = [];
blana.eventValues = {{'rlCL','CL'}};
blana.eventValuesSplit = {{ {'rlCLrtcrt'} {'rFix'} }};
blana.trl_expr = {{ 
                    {'sttp==2 & crct==1'} ... %{'rlCLrtcrt'}  
                    {'sttp==1 | sttp==2'}; 
                 }};        
blana.trl_order = ana.trl_order;             
inana{2} = blana;

blana = [];
blana.eventValues = {{'rlOL','OL'}};
blana.eventValuesSplit = {{ {'rlOLrt'} {'rFix'} }};
blana.trl_expr = {{ {'sttp==1'} ... %{'rlOLrtcrt'}
                    {'sttp==1 | sttp==2'}; 
                 }};     
blana.trl_order = ana.trl_order;
inana{3} = blana;

blana = [];
blana.eventValues = {{'rlCL','CL'}};
blana.eventValuesSplit = {{ {'rlCLrt'} {'rFix'} }};
blana.trl_expr = {{ 
                    {'sttp==2'} ... %{'rlCLrtcrt'}  
                    {'sttp==1 | sttp==2'}; 
                 }};        
blana.trl_order = ana.trl_order;
inana{4} = blana;

blana = [];
blana.eventValues = {{'OL','OL'}};
blana.eventValuesSplit = {{ {'OLrt'} {'rFix'} }};
blana.trl_expr = {{ 
                    {'sttp==1'} ... %{'rlCLrtcrt'}  
                    {'sttp==1 | sttp==2'}; 
                 }};        
blana.trl_order = ana.trl_order;
inana{5} = blana;

blana = [];
blana.eventValues = {{'CL','CL'}};
blana.eventValuesSplit = {{ {'CLrt'} {'rFix'} }};
blana.trl_expr = {{ 
                    {'sttp==2'} ... %{'rlCLrtcrt'}  
                    {'sttp==1 | sttp==2'}; 
                 }};        
blana.trl_order = ana.trl_order;
inana{6} = blana;

blana = [];
blana.eventValues = {{'OL','OL'}};
blana.eventValuesSplit = {{ {'OLrtcrt'} {'rFix'} }};
blana.trl_expr = {{ 
                    {'sttp==1 & crct==1'} ... %{'rlCLrtcrt'}  
                    {'sttp==1 | sttp==2'}; 
                 }};        
blana.trl_order = ana.trl_order;
inana{7} = blana;

blana = [];
blana.eventValues = {{'CL','CL'}};
blana.eventValuesSplit = {{ {'CLrtcrt'} {'rFix'} }};
blana.trl_expr = {{ 
                    {'sttp==2 & crct==1'} ... %{'rlCLrtcrt'}  
                    {'sttp==1 | sttp==2'}; 
                 }};        
blana.trl_order = ana.trl_order;
inana{8} = blana;

%high vs. low conf
blana = [];
blana.eventValues = {{'rlOL','rlCL','OL','CL'}};
blana.eventValuesSplit = {{ {'rlOLrthcf'} {'rlCLrthcf'} {'OL'} {'CL'} }};
blana.trl_expr = {{ 
                    {'sttp==1 & conf>2 & crct==1'} ... 
                    {'sttp==2 & conf>2 & crct==1'} ...
                    {'sttp==1'}...
                    {'sttp==2'}...
                 }};        
blana.eventValuesComb = {{ {'rlOLrthcf', 'rlCLrthcf'} {'OL','CL'}}};
blana.eventValuesNew = {{ 'rlRetrthcfcrct', 'rFix' }};             
blana.trl_order = ana.trl_order;
inana{9} = blana;

blana = [];
blana.eventValues = {{'rlOL','rlCL','OL','CL'}};
blana.eventValuesSplit = {{ {'rlOLrtlcf'} {'rlCLrtlcf'} {'OL'} {'CL'}}};
blana.trl_expr = {{ 
                    {'sttp==1 & conf<3 & crct==1'} ... 
                    {'sttp==2 & conf<3 & crct==1'} ...
                    {'sttp==1'}...
                    {'sttp==2'}...
                 }};        
blana.eventValuesComb = {{ {'rlOLrtlcf', 'rlCLrtlcf'} {'OL','CL'} }};
blana.eventValuesNew = {{ 'rlRetrtlcfcrct', 'rFix' }};
blana.trl_order = ana.trl_order;
inana{10} = blana;


% rest of defs

cfg_bl = [];

trgCond = {'rlOLrtcrt','rlCLrtcrt','rlOLrt','rlCLrt','OLrt','CLrt','OLrtcrt','CLrtcrt','rlRetrthcfcrct','rlRetrtlcfcrct'};
blCond = repmat({'rFix'},size(trgCond));

cfg_bl.loadcfg = deevLoadConfig;
cfg_bl.loadcfg.ftype = 'fourier';
cfg_bl.loadcfg.loadMethod = 'trialinfo';
cfg_bl.loadcfg.keeptrials = 'yes';
cfg_bl.loadcfg.baseline_type = '';
cfg_bl.loadcfg.norm_trials = 'average'; %watch out for 'norm_trial' gotcha
cfg_bl.loadcfg.baseline_time = [];

cfg_bl.keeptrials = 'no';
cfg_bl.method = 'zscore';
cfg_bl.twin = [-.3 -.1];
cfg_bl.param = 'powspctrm';

cfg_bl.nullRT = true;
cfg_bl.bufRT = .3;

for icond = 9:10%length(trgCond)
    cfg_bl.colinds = true(size(ana.trl_order.rlCL));%they're all the same size anyway
    cfg_bl.colinds(1) = 0; %condn col can't be used in matching
    cfg_bl.blCond = blCond{icond};
    cfg_bl.trgCond = trgCond{icond};
    
    [bldata, exper] = deevTrialBL(exper,dirs,inana{icond},cfg_bl);
    data_pow.ses1.(['tbl' cfg_bl.trgCond]) = bldata;
end
%add new cond to eventValues
ana.eventValues = {fieldnames(data_pow.ses1)};
%select bl corrected data
ind = ~cellfun(@isempty,regexp(ana.eventValues{1},'^tbl.*'));
ana.eventValues = {ana.eventValues{1}(ind)};
   

%% find low trial counts
%cond baseline
conds = fieldnames(exper.nTrials.ses1);
nTrials = [];
for icond = 1:length(conds)
    nTrials = cat(2,nTrials,exper.nTrials.ses1.(conds{icond}));
end
%%
%trial baseline
conds = fieldnames(exper.nTrials);
%conds = {'OL','CL'};
nTrials = [];
for icond = 1:length(conds)
    nTrials = cat(2,nTrials,exper.nTrials.(conds{icond}));
end
%%
nthres = 20;
exper.badSub = sum(nTrials<nthres,2)>1;

ana = mm_ft_elecGroups(ana);


%% find bad performing subjects
cfg_beh = [];
cfg_beh.dir = fullfile(dirs.dataroot,dirs.behDir);
cfg_beh.filefiltstr = 'deevEEGExp_1EEG_4blocks_36events_Sub[0-9]{1,3}\.txt';
cfg_beh.badsubs = [100];
out = deevBlkGetLogDep(cfg_beh);

ana = mm_ft_elecGroups(ana);

%% get the grand average

% set up strings to put in grand average function
cfg_ana = [];
cfg_ana.is_ga = 0;
cfg_ana.data_str = 'data_pow';

cfg_ft = [];
%cfg_ft.keepindividual = 'no';
cfg_ft.keepindividual = 'yes';

if strcmp(cfg_ana.data_str,'data_pow')
  ga_pow = struct;
elseif strcmp(cfg_ana.data_str,'data_pow_log')
  ga_pow_log = struct;
elseif strcmp(cfg_ana.data_str,'data_coh')
  ga_coh = struct;
elseif strcmp(cfg_ana.data_str,'data_evoked')
  ga_evoked = struct;
end

for ses = 1:length(exper.sesStr)

    cfg_ana.conditions = ana.eventValues{ses};
    cfg_ana.sub_str = mm_catSubStr_multiSes(cfg_ana,exper,ses);
    
    for evVal = 1:length(ana.eventValues{ses})
      %tic
      fprintf('Running ft_freqgrandaverage on %s...',ana.eventValues{ses}{evVal});
      if strcmp(cfg_ana.data_str,'data_pow')
        cfg_ft.parameter = 'powspctrm';
        ga_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{ses}{evVal})));
      elseif strcmp(cfg_ana.data_str,'data_pow_log')
        cfg_ft.parameter = 'powspctrm';
        ga_pow_log.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{ses}{evVal})));
      elseif strcmp(cfg_ana.data_str,'data_coh')
        %cfg_ft.parameter = 'plvspctrm';
        cfg_ft.parameter = 'powspctrm';
        ga_coh.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{ses}{evVal})));
      elseif strcmp(cfg_ana.data_str,'data_evoked')
        cfg_ft.parameter = 'powspctrm';
        ga_evoked.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{ses}{evVal})));
      end
      fprintf('Done.\n');
      %toc
  end
end


%% plot power 
cfg_ft = [];
%cfg_ft.xlim = [0.5 6];
%cfg_ft.ylim = [5.8 6.2];
%cfg_ft.ylim = [9.8 10.2];
%cfg_ft.ylim = [19.8 20.2]; 
%cfg_ft.ylim = [3 30];
%cfg_ft.ylim = [3 8];
%cfg_ft.ylim = [8 12];
%cfg_ft.ylim = [12 28];
%cfg_ft.ylim = [28 64];
%cfg_ft.zlim = [-1 1];
cfg_ft.zlim = [-1.5 1.5];
%elseif strcmp(cfg_ft.baselinetype,'relative')
%cfg_ft.zlim = 'maxmin';
%end
cfg_ft.showlabels = 'yes';
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.layout = ft_prepare_layout([],ana);
cfg_ft.roi = ana.elecGroups{ismember(ana.elecGroupsStr,'noEyeABH')};

for ses = 1:length(ana.eventValues)
    for evVal = 1:length(ana.eventValues{ses})
      figure
      ft_singleplotTFR(cfg_ft,ga_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal})); 
      hold on
%       entn = regexp(ana.eventValues{ses}{evVal},'([0-9]+)$','tokens');
%       entn = str2num(entn{1}{1});
%       if entn>0
%           plot(xlim,[entn entn],'--k','linewidth',2)
%       end
      %ft_topoplotTFR(cfg_ft,ga_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}));              
      %ft_multiplotTFR(cfg_ft,ga_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}));
      set(gcf,'Name',sprintf('%s',strrep(ana.eventValues{ses}{evVal},'_','-')))
    end

end

%% plot the contrasts


cfg_plot = [];
cfg_plot.plotTitle = 1;
figure('color','white');
% comparisons to make
cfg_ft = [];
cfg_ft.parameter = 'powspctrm';
cfg_ft.zlim = [-1 1]; % pow
%cfg_ft.zlim = 'maxmin'; % pow
cfg_ft.interactive = 'yes';
cfg_ft.colormap = parula;
cfg_ft.colorbar = 'no';
cfg_plot.zlabel = 'relchange power';
%%%%%%%%%%%%%%%
% Type of plot
%%%%%%%%%%%%%%%

%cfg_plot.conditions = {{'tblrlRetrthcf','tblrlRetrtlcf'}};
cfg_plot.conditions = {{'varOL','varCL'}};
%cfg_plot.conditions = {{'tblrlCLrtcrt','tblrlOLrtcrt'}};
%cfg_plot.conditions = {{'CL','OL'}};
cfg_plot.ftFxn = 'ft_singleplotTFR';
cfg_plot.roi = {'noEyeABH'};
%cfg_plot.ftFxn = 'ft_topoplotTFR';
cfg_ft.marker = 'off';
cfg_ft.comment = 'no';
cfg_plot.subplot = 0;
cfg_plot.plotTitel = 0;
%cfg_ft.xlim = [-2.5 -2]; % time
cfg_ft.ylim = [3 50]; %freq
%cfg_ft.avgovertime = 'yes';
%cfg_ft.xlim = (-3:0.5:0); % time
%cfg_ft.xlim = [-.3 3]; % time


%cfg_plot.ftFxn = 'ft_multiplotTFR';
cfg_ft.showlabels = 'no';
% cfg_ft.comment = '';

mm_ft_contrastTFR(cfg_ft,cfg_plot,ana,exper,files,dirs,ga_pow,1);

%% cluster stat




%% line plot

files.saveFigs = false;
files.figPrintFormat = 'png';

cfg = [];
cfg.parameter = 'powspctrm';
t = [0 3]; intv = .04;
cfg.times = [t(1):intv:t(2)-intv; t(1)+intv:intv:t(2)]';
%cfg.times = [-0.5:0.02:5.98; -0.48:0.02:6.0]';
%cfg.times = [.90:0.02:1.30; .92:0.02:1.32]';
%cfg.times = [0 1];

cfg.conditions = {{'tblCLrtcrt','tblOLrtcrt'}};
%cfg.conditions = {{'tblrlCLrtcrt','tblrlOLrtcrt'}};

cfg.graphcolor = 'gbgbcc';
cfg.linestyle = {'-','-','--','--'};

f = mm_freqSet();
% cfg.freqs = ana.freq.theta;
% cfg.freqs = ana.freq.alpha_lower;
% cfg.freqs = ana.freq.alpha_upper;
% cfg.freqs = ana.freq.beta_lower;
%cfg_ana.freqs = {[5.8 6.2],[9.8 10.2],[19.8 20.2],[39.8 40.2]};
%cfg_ana.freqs = {[3 8],[8 12],[12 30]};
cfg_ana.freqs = {f.theta,f.alpha,f.beta};
%cfg.rois = {sigElecs};

cfg.plotTitle = false;
cfg.plotLegend = true;
cfg.legendloc = 'NorthEast';

cfg.plotErrorBars = true;
cfg.eb_transp = true;

%cfg.plotClusSig = true;
% cfg.clusAlpha = 0.1;
%cfg.clusTimes = cfg_ana.latencies;
% % cfg.clusTimes = [-0.2:0.2:0.8; 0:0.2:1.0]';
% %cfg.clusTimes = [-0.18:0.1:0.92; -0.1:0.1:1.0]'; % 100 no overlap
% cfg.clusTimes = [-0.18:0.2:0.92; 0:0.2:1.0]'; % 200 no overlap
% cfg.clusLimits = false;

cfg.linewidth = 2;
% cfg.limitlinewidth = 0.5;
% cfg.textFontSize = 10;

%cfg.yminmax = [-2 2];
%cfg.yminmax = [-0.6 0.6];
%cfg.yminmax = [-0.5 0.2];
cfg.nCol = 3;

cfg.rois = {'LAS','PS'};

% % whole power
% cfg.type = 'line_pow';
% cfg.clusDirStr = '_zpow_-400_-200';
% cfg.ylabel = 'Z-Trans Pow';
for ifreq = 1:length(cfg_ana.freqs)
    cfg.freqs = cfg_ana.freqs{ifreq};
    mm_ft_lineTFR(cfg,ana,exper,files,dirs,ga_pow);
end

%% plot avg pow over time window 

cfg.latency = [-3 -1];
cfg.avgovertime = 'yes';

%cfg_ana.frequency = {[5.8 6.2],[9.8 10.2],[19.8 20.2]};
cfg_ana.frequency = {[3 8],[8 12],[8 20]};
cfg.avgoverfreq = 'yes';

roi = 'PS';

roiind = ismember(ana.elecGroupsStr,roi);
cfg.channel = ana.elecGroups{roiind};
cfg.avgoverchan = 'yes';
cfg.nanmean = 'yes';
avgdata = [];
cbardata = {};
for ifreq = 1:length(cfg_ana.frequency)
    
    cfg.frequency = cfg_ana.frequency{ifreq};
    conds = {'tblrlCLrtcrt','tblrlOLrtcrt'};
    %conds = {'CL','OL'};
    %get avg data for each condition
    bardata = [];
    for icond = 1:length(conds)
        for isub = 1:length(exper.subjects)
            avgdata(isub).(conds{icond}) =  ft_selectdata(cfg, data_pow.ses1.(conds{icond}).sub(isub).data);
            bardata(isub,icond) = avgdata(isub).(conds{icond}).powspctrm;
        end
    end
    figure('color','white');
    hold on    
    bar(mean(bardata));
    plot(bardata','.-','markersize',30);
    box off
    set(gca,'Xtick',1:4);
    set(gca,'Xticklabel',conds);
    ylabel(sprintf('%.01f-%.01f avg power',cfg.frequency(1), cfg.frequency(2)));
    set(gca,'fontsize',24);
    cbardata{ifreq} = bardata;
end
%%
%correlation of memory with entrained power
icond = 2;
dpcdval = nan(length(conds)-1,length(cfg_ana.frequency));
dpcval = nan(length(conds),length(cfg_ana.frequency));

pcdval = nan(length(conds)-1,length(cfg_ana.frequency));
pcval = nan(length(conds),length(cfg_ana.frequency));
for ifreq = 1:length(cfg_ana.frequency) 
    
    pw = cbardata{ifreq};
    dpw = nan(size(pw,1),size(pw,2)-1);
    for icol = 2:size(pw,2)
        dpw(:,icol-1) = pw(:,icol)-pw(:,1);
    end
    dp = out.results.bycond.dprime{icond}';
    ph = out.results.bycond.pairhit{icond}';
    
    dpcdval(:,ifreq) = corr(dpw,dp);
    dpcval(:,ifreq) = corr(pw,dp);
    
    pcdval(:,ifreq) = corr(dpw,ph);
    pcval(:,ifreq) = corr(pw,ph);
    
end



