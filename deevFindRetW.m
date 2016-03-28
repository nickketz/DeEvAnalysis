function out = deevFindRetW(substr, condN, chansel, adfile)
% function to estimate optimal time x freq window for detecting retrieval
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

%% config for loading data
cfg = deevLoadConfig;

%% do simulated annealing search for condN eventValue

ana = mm_ft_elecGroups(ana);
chanstr = ana.elecGroups{strcmp(ana.elecGroupsStr,chansel)};

%load training data
% load OL and CL trials
ana.eventValues = {{'OL','CL'}};

%load data and do preprocessing
data_pow =  mm_ft_loadData_multiSes(cfg,exper,dirs,ana);
cfg_app = []; cfg_app.parameter = 'powspctrm'; cfg_app.appenddim = 'rpt';
data = ft_appendfreq(cfg_app,data_pow.ses1.OL.sub.data,data_pow.ses1.CL.sub.data);
trialinfo = data.trialinfo;
clear data_pow;

%
%other needed vars
chaninds = ismember(data.label,chanstr);
assignin('base','chaninds',chaninds);

%lower and upper bounds of search
lb = [0 .1 data.freq(1) 1];%lower bounds

t0 = find(data.time==0);
tnan = squeeze(isnan(data.powspctrm(1,chaninds,1,:)));
tnan = sum(tnan,1)>0;
tnan(1:t0) = 0;
[~,~,tmax] = findnearest(min(data.time(tnan))-lb(2), data.time,-1);
[~,~,fmax] = findnearest(data.freq(end)-lb(4), data.freq,-1);
%upper bounds
ub = [  tmax ... %stay away from boundary nans 
        2 ... %max of 2sec time window
        fmax - lb(4)... %use whole window minus lowerbound of freq window
        data.freq(end)/2]; %max window 1/2 of total
    

%trim dpos and assign in base
trlord = ana.trl_order.OL; 
if (size(data.trialinfo,2)-length(trlord)==1) & (trialinfo(1,2)==str2double(exper.subjects{1}(end-1:end)))
   trlord = cat(2,'condn',trlord); %adding condn onto trlord 
end
ordncol = strcmp('ordn',trlord);
sttp = trialinfo(:,ordncol);
cuetp = (sttp-mod(sttp,10))/10;
targtp = mod(sttp,10);

%select pos and neg trials
cfg_sel = [];
cfg_sel.trials = cuetp==condN;
[~,~,tmax] = findnearest(min(data.time(tnan)),data.time);
cfg_sel.latency = [0 tmax];
postrialinfo = trialinfo(cfg_sel.trials,:);
assignin('base','dpos',ft_selectdata(cfg_sel,data));

cfg_sel.trials = cuetp~=condN;
assignin('base','dneg',ft_selectdata(cfg_sel,data));

%assign other needed stuff in base
assignin('base','ub',ub);
assignin('base','substr',substr);
assignin('base','evtStr',conds{condN});
assignin('base','slide',slide);

%search params
fun = @deev_tfWObjective;
x0 = [0 1 4 10];%startTime, timeWindow, lowerFreq, freqWindow

options = saoptimset('PlotFcns',{@deev_saplotbestx, @saplotbestf, @saplotf, @saplottemperature},...
   'TemperatureFcn',@temperatureboltz);%'ReannealInterval',50, 'InitialTemperature',[500 500 500 500]);

[x,fmin,ef,output] = simulannealbnd(fun,x0,lb,ub,options);
%x = x0; fmin = 1; ef = nan; output = nan;
%get optimal template based on search params
[W,t,f] = deevGetW(x);


%% cross-validate W
[avgtestsim,avgtrainsim,twin] = deevCrosValW(x,slide,evalin('base','dpos'),-1);

%clear large data sets from base ws
evalin('base','clear dneg dpos');

%% eval W on negative condition
% cfg_sel.trials = cuetp~=condN & targtp~=condN;
% dneg = ft_selectdata(cfg_sel,data);
dneg = evalin('base','dneg');
negtrialinfo = dneg.trialinfo;
%get similarity score
[negsim, negtsamps] = deevEvalW(W,t,f,chaninds,dneg);


%% save
% name: W_condN_chansel_slide.mat in
% dir:  dirs.saveDirProc/substr/ses1/
savename = sprintf('WRet_%s_%s_%0.2f.mat',conds{condN},chansel,slide);
savename = fullfile(dirs.saveDirProc,substr,'ses1',savename);

savevars = {'substr','condN','conds','chansel','slide'...
    'W','t','f','x','fmin','ef','output','cfg',...
    'avgtestsim','avgtrainsim','twin','postrialinfo',...
    'negsim','negtsamps','negtrialinfo'};

save(savename,savevars{:});

