%preliminary power analysis to check on flckr condtions

%load analysis details
%adFile = '/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/v1/ft_data/OL_CL_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
%adFile = '/home/nike3851/work/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/v1/ft_data/OL_CL_encLoc_encPer_encObj_encAni_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
adFile = '/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/v1/ft_data/OL_CL_encLoc_encPer_encObj_encAni_encNotLoc_encNotPer_encNotObj_encNotAni_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
%adFile = '/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/v1/ft_data/OL_CL_encLoc_encPer_encObj_encAni_encNotLoc_encNotPer_encNotObj_encNotAni_eq0_art_ftManual/pow_wavelet_w4_pow_3_50/analysisDetails.mat';

load(adFile)

%%
%remove subs?
conds = exper.eventValues;
tmpexper = exper;
%find subs with 0 trials in any condition
nTrials = [];
for icond = 1:length(conds)
    nTrials = cat(2,nTrials,exper.nTrials.(conds{icond}));
end
%badInd = sum(nTrials<20,2);
%also remove sub20 for no behav data
badInd(~strcmp(exper.subjects,'DeEv_sub11')) = 1;

exper = nk_rmSubs(exper,badInd);

%%
%define trials of interest
ana.eventValues = {{'OL','CL'}};
ana.eventValues = {exper.eventValues};

%load data into ft data struct
%[data_pow,exper] = mm_loadSubjectData(exper,dirs,ana,'pow',keeptrials,'trialinfo',false);
cfg = [];

cfg.loadMethod = 'seg';
%cfg.loadMethod = 'trialinfo';
cfg.latency = 'all';
cfg.frequency = 'all';

cfg.keeptrials = 'no';
%cfg.keeptrials = 'yes';
cfg.equatetrials = 'no';
% %cfg.equatetrials = 'yes';

cfg.rmPreviousCfg = true;

% type of input (used in the filename to load)
cfg.ftype = 'pow';
% cfg.ftype = 'fourier';

% type of output: 'pow', 'coh', 'phase'
cfg.output = 'pow';

% transformation: 'log10', 'log', 'vec'
cfg.transform = '';
% cfg.transform = 'log10';
% cfg.transform = 'vec';

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


[data_pow,exper] =  mm_ft_loadData_multiSes(cfg,exper,dirs,ana);

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
cfg_ft.keepindividual = 'no';
%cfg_ft.keepindividual = 'yes';

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
cfg_ft.xlim = [-1 2];
%cfg_ft.ylim = [5.8 6.2];
%cfg_ft.ylim = [9.8 10.2];
%cfg_ft.ylim = [19.8 20.2];
%cfg_ft.ylim = [3 30];
%cfg_ft.ylim = [3 8];
%cfg_ft.ylim = [8 12];
%cfg_ft.ylim = [12 28];
%cfg_ft.ylim = [28 64];
%cfg_ft.zlim = [-1 1];
%cfg_ft.zlim = [-2 2];
%elseif strcmp(cfg_ft.baselinetype,'relative')
cfg_ft.zlim = [-.5 .5];
%end
cfg_ft.showlabels = 'yes';
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.layout = ft_prepare_layout([],ana);
%cfg_ft.channel = ana.elecGroups{ismember(ana.elecGroupsStr,'PS2')};

for ses = 1:length(ana.eventValues)

    for evVal = 1:length(ana.eventValues{ses})
      figure
      ft_singleplotTFR(cfg_ft,ga_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal})); 
      hold on
      entn = regexp(ana.eventValues{ses}{evVal},'([0-9]+)$','tokens');
      entn = str2num(entn{1}{1});
      if entn>0
          plot(xlim,[entn entn],'--k','linewidth',2)
      end
      %ft_topoplotTFR(cfg_ft,ga_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}));              
      %ft_multiplotTFR(cfg_ft,ga_pow.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}));
      set(gcf,'Name',sprintf('%s',strrep(ana.eventValues{ses}{evVal},'_','-')))
    end

end

%% plot the contrasts


cfg_plot = [];
cfg_plot.plotTitle = 1;

% comparisons to make
cfg_plot.conditions = {{'flckr6','flckr0'},{'flckr10','flckr0'},{'flckr20','flckr0'}};
%cfg_plot.conditions = {{'flckr6','flckr0'},{'flckr6','flckr10'},{'flckr6','flckr20'}};
%cfg_plot.conditions = {{'flckr20','flckr0'}};
cfg_ft = [];
%cfg_ft.xlim = [-1 2]; % time
%cfg_ft.ylim = [4 8]; % freq
%cfg_ft.ylim = [8 12]; % freq
%cfg_ft.ylim = [9.8 10.2];
%cfg_ft.ylim = [5.8 6.2]; % freq
cfg_ft.ylim = [19.8 20.2]; % freq
%cfg_ft.ylim = [3 30]; % freq
%cfg_ft.ylim = [28 50]; % freq
cfg_ft.parameter = 'powspctrm';
cfg_ft.zlim = [-.5 .5]; % pow

%cfg_ft.colorbar = 'no';
cfg_ft.interactive = 'yes';
cfg_ft.colormap = parula;
cfg_ft.colorbar = 'no';
cfg_plot.zlabel = 'relchange power';
%%%%%%%%%%%%%%%
% Type of plot
%%%%%%%%%%%%%%%

%cfg_plot.ftFxn = 'ft_singleplotTFR';
cfg_plot.ftFxn = 'ft_topoplotTFR';
%cfg_ft.marker = 'on';
cfg_ft.marker = 'labels';
cfg_ft.markerfontsize = 9;
cfg_ft.comment = 'no';
%cfg_ft.xlim = [0.5 0.8]; % time
cfg_plot.subplot = 0;
cfg_ft.xlim = [0 1]; % time
cfg_ft.avgovertime = 'no';
%cfg_ft.xlim = (0:0.05:1.0); % time
cfg_plot.roi = {'PS2'};


% cfg_plot.ftFxn = 'ft_multiplotTFR';
cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';

mm_ft_contrastTFR(cfg_ft,cfg_plot,ana,exper,files,dirs,ga_pow,1);


%% line plot

files.saveFigs = false;
files.figPrintFormat = 'png';

cfg = [];
cfg.parameter = 'powspctrm';

cfg.times = [-0.5:0.02:1.98; -0.48:0.02:2.0]';
%cfg.times = [.90:0.02:1.30; .92:0.02:1.32]';
%cfg.times = [0 1];

%cfg.conditions = {{'flckr6','flckr0'},{'flckr10','flckr0'},{'flckr20','flckr0'}};
cfg.conditions = {{'flckr0','flckr6','flckr10','flckr20'}};
% cfg.rename_conditions = {{'Space12 P2 Recalled','Space12 P2 Forgot','Space32 P2 Recalled','Space32 P2 Forgot'}};

cfg.graphcolor = 'bcrm';
cfg.linestyle = {'-','--','-','--'};

% cfg.freqs = ana.freq.theta;
% cfg.freqs = ana.freq.alpha_lower;
% cfg.freqs = ana.freq.alpha_upper;
% cfg.freqs = ana.freq.beta_lower;
cfg_ana.freqs = {[5.8 6.2],[9.8 10.2],[19.8 20.2],[39.8 40.2]};
%cfg_ana.freqs = {[4 8],[8 12],[12 30]};
% cfg.rois = {sigElecs};

cfg.plotTitle = false;
cfg.plotLegend = true;
cfg.legendloc = 'NorthEast';

cfg.plotErrorBars = false;
cfg.eb_transp = true;

cfg.plotClusSig = false;
% cfg.clusAlpha = 0.1;
% %cfg.clusTimes = cfg.times;
% % cfg.clusTimes = [-0.2:0.2:0.8; 0:0.2:1.0]';
% %cfg.clusTimes = [-0.18:0.1:0.92; -0.1:0.1:1.0]'; % 100 no overlap
% cfg.clusTimes = [-0.18:0.2:0.92; 0:0.2:1.0]'; % 200 no overlap
% cfg.clusLimits = false;

cfg.linewidth = 2;
% cfg.limitlinewidth = 0.5;
% cfg.textFontSize = 10;

cfg.yminmax = [-2 2];
%cfg.yminmax = [-0.6 0.6];
%cfg.yminmax = [-0.5 0.2];
cfg.nCol = 3;

cfg.rois = 'PS2';

% % whole power
% cfg.type = 'line_pow';
% cfg.clusDirStr = '_zpow_-400_-200';
% cfg.ylabel = 'Z-Trans Pow';
for ifreq = 1:length(cfg_ana.freqs)
    cfg.freqs = cfg_ana.freqs{ifreq};
    mm_ft_lineTFR(cfg,ana,exper,files,dirs,ga_pow);
end

%% plot avg pow over time window 

cfg.latency = [0.2 1];
cfg.avgovertime = 'yes';

cfg_ana.frequency = {[5.8 6.2],[9.8 10.2],[19.8 20.2]};
%cfg_ana.frequency = {[3 8],[8 12],[12 30]};
cfg.avgoverfreq = 'yes';

roi = 'PS2';

roiind = ismember(ana.elecGroupsStr,roi);
cfg.channel = ana.elecGroups{roiind};
cfg.avgoverchan = 'yes';
avgdata = [];
cbardata = {};
for ifreq = 1:length(cfg_ana.frequency)
    
    cfg.frequency = cfg_ana.frequency{ifreq};
    conds = {'flckr0','flckr6','flckr10','flckr20'};
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
    plot(bardata','.','markersize',30);
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



