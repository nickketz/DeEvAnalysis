function out = deevFindW(substr, condN, chansel, adfile)
% function to estimate optimal time x freq window for detecting encoding
% stim categories
% can be run through sbatch where each array is specified to condN
%   i.e. sbatch -a 1:4 runFindW.sh 
%
%
% input:
%   substr: string of subject to run 
%   condN: int, condition number where conditions= {'encLoc','encPer','encObj','encAni'};
%   chansel: str, which elecGroupStr to use to sub-select chans in search
%   slide: temporal slide between evaluations in seconds, def .1
%   adfile: str, analysis details file
%
% output
%   saves output file W_condN_chansel.mat in
%   dirs.saveDirProc/substr/ses1/
%       includes W, x, fmin, exitflag, output log of search, and eval of W
%       in reference to OL and CL conditions as well as the negative
%       encoding condition (i.e. searched encPer, eval on encNotPer)



%old output
%   out: output structure with fields
%       s: array, optimal start time for each eventValue
%       e: array of optimal end times for each eventValue
%       l: array, optimal lower freqs x eventValue
%       u: array, optimal upper freqs x eventValue
%       fmin: array, min objective achieved for each eventValue
%

%TODO add cfg as input with defaults

%analysis details file
if ~exist('adfile','var') 
    adFile = '/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/v1/ft_data/OL_CL_encLoc_encPer_encObj_encAni_encNotLoc_encNotPer_encNotObj_encNotAni_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
%    adFile = '/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/long/ft_data/OL_CL_encLoc_encPer_encObj_encAni_encNotLoc_encNotPer_encNotObj_encNotAni_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
end

% elecGroupStr to sub-select chans for search
if ~exist('chansel','var')
    chansel = 'noEyeABH';
end

if ~exist('slide','var')
    slide = .1;
end

load(adFile) %contains exper, ana, dir structures

%select specific subject
badInd = zeros(size(exper.subjects));
badInd(~strcmp(exper.subjects,substr)) = 1;
exper = nk_rmSubs(exper,badInd);

% condition names to assign per condN input
conds = {'encLoc','encPer','encObj','encAni'};

%% load data
cfg = [];

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

%% do simulated annealing search for each eventValue

ana = mm_ft_elecGroups(ana);
chanstr = ana.elecGroups{strcmp(ana.elecGroupsStr,chansel)};

%search condition
evtStr = conds{condN};
ana.eventValues = {{evtStr}};
%load data and do preprocessing
data_pow =  mm_ft_loadData_multiSes(cfg,exper,dirs,ana);

%other needed vars
chaninds = ismember(data_pow.ses1.(evtStr).sub.data.label,chanstr);
assignin('base','chaninds',chaninds);

%lower and upper bounds of search
lb = [0 .1 data_pow.ses1.(evtStr).sub.data.freq(1) 1];%lower bounds

t0 = find(data_pow.ses1.(evtStr).sub.data.time==0);
tnan = squeeze(isnan(data_pow.ses1.(evtStr).sub.data.powspctrm(1,chaninds,1,:)));
tnan = sum(tnan,1)>0;
tnan(1:t0) = 0;
[~,~,tmax] = findnearest(min(data_pow.ses1.(evtStr).sub.data.time(tnan))-lb(2), data_pow.ses1.(evtStr).sub.data.time,-1);
[~,~,fmax] = findnearest(data_pow.ses1.(evtStr).sub.data.freq(end)-lb(4), data_pow.ses1.(evtStr).sub.data.freq,-1);
%upper bounds
ub = [  tmax ... %stay away from boundary nans 
        2 ... %max of 2sec time window
        fmax - lb(4)... %use whole window minus lowerbound of freq window
        data_pow.ses1.(evtStr).sub.data.freq(end)/2]; %max window 1/2 of total
    
%trim dpos and assign in base
[~,~,tmax] = findnearest(min(data_pow.ses1.(evtStr).sub.data.time(tnan)),data_pow.ses1.(evtStr).sub.data.time);
cfg_sel = [];
cfg_sel.latency = [0 tmax];
assignin('base','dpos',ft_selectdata(cfg_sel,data_pow.ses1.(evtStr).sub.data));


% load negative condition
negevtStr = [evtStr(1:3) 'Not' evtStr(4:end)];
ana.eventValues = {{negevtStr}};
%load data and do preprocessing
data_pow =  mm_ft_loadData_multiSes(cfg,exper,dirs,ana);
%subselect first two presentations or OL conditions to avoid latent
%reactivation of non-target stims
trlord = ana.trl_order.(negevtStr);
if size(trlord,2) ~= size(data_pow.ses1.(negevtStr).sub.data.trialinfo,2)
    trlord = cat(2,'condn',trlord);%add condn as first col of trlord
end
colind = strcmp('blkn',trlord);
cfg.trials =  data_pow.ses1.(negevtStr).sub.data.trialinfo(:,colind)~=3;
if sum(cfg.trials)<3
    warning('too few trials, using all blocks instead of all but blk 3');
    cfg.trials = ones(size(cfg.trials));
end
assignin('base','dneg',ft_selectdata(cfg_sel, data_pow.ses1.(negevtStr).sub.data));

%assign other needed stuff in base
assignin('base','ub',ub);
assignin('base','substr',substr);
assignin('base','evtStr',evtStr);
assignin('base','slide',slide);

%search params
fun = @deev_tfWObjective;
x0 = [0 1 4 10];%startTime, timeWindow, lowerFreq, freqWindow

options = saoptimset('PlotFcns',{@deev_saplotbestx, @saplotbestf, @saplotf, @saplottemperature},...
   'TemperatureFcn',@temperatureboltz);%'ReannealInterval',50, 'InitialTemperature',[500 500 500 500]);

[x,fmin,ef,output] = simulannealbnd(fun,x0,lb,ub,options);

%get optimal template based on search params
[W,t,f] = deevGetW(x);

%clear large data sets from base ws
evalin('base','clear dneg dpos');


%% eval W on negative condition

%get similarity score
[negsim, negtsamps] = deevEvalW(W,t,f,chaninds,data_pow.ses1.(negevtStr).sub.data);


%% eval on OL trials
ana.eventValues = {{'OL'}};

%load data and do preprocessing
data_pow =  mm_ft_loadData_multiSes(cfg,exper,dirs,ana);

%get similarity score
[olsim, oltsamps] = deevEvalW(W,t,f,chaninds,data_pow.ses1.OL.sub.data);
oltrialinfo = data_pow.ses1.OL.sub.data.trialinfo;

%% eval on CL
ana.eventValues = {{'CL'}};

%load data and do preprocessing
data_pow =  mm_ft_loadData_multiSes(cfg,exper,dirs,ana);

%get similarity score
[clsim, cltsamps] = deevEvalW(W,t,f,chaninds,data_pow.ses1.CL.sub.data);
cltrialinfo = data_pow.ses1.CL.sub.data.trialinfo;

%% save
% name: W_condN_chansel_slide.mat in
% dir:  dirs.saveDirProc/substr/ses1/
clear data_pow
savename = sprintf('W_%s_%s_%0.2f.mat',conds{condN},chansel,slide);
savename = fullfile(dirs.saveDirProc,substr,'ses1',savename);

savevars = {'substr','condN','conds','chansel','slide'...
    'W','t','f','x','fmin','ef','output','cfg',...
    'negsim','negtsamps','olsim','oltsamps','clsim','cltsamps','oltrialinfo','cltrialinfo'};

save(savename,savevars{:});

