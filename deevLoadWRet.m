function [wfile,subjplot] = deevLoadWRet(substr,cfg)
%
% function to get weight vector W for a given subject train on retrieval 
% conditions, 
%
% input:
%   substr: string subject to retrieve W info
%   cfg: config struct with optional fields
%       doplots: bool, plot histograms of performance
%       adFile: string, analysis details file, def:'/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/v1/ft_data/OL_CL_encLoc_encPer_encObj_encAni_encNotLoc_encNotPer_encNotObj_encNotAni_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat';
%       chanstr: string, elecGroupStr used to select W, def: 'noEyeABH'
%       slide: float, slide used to eval W, def: 0.1
%       conds: cell, conditions to eval W on, def:  {'encLoc','encPer','encObj','encAni'}
%       subplot: handel into shared plot across subjects
%
% output:
%   wfile: struct of data from W file saved from call to deevFindW
%   subplot: handel into figure for shared plotting
%

if ~exist('cfg','var'),             cfg = [];                           end
if ~isfield(cfg,'doplots'),         cfg.doplots = 0;                    end
if ~isfield(cfg,'adFile'),          cfg.adFile = '/work/ccnlab/users/nike3851/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/v1/ft_data/OL_CL_encLoc_encPer_encObj_encAni_encNotLoc_encNotPer_encNotObj_encNotAni_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat'; end
if ~isfield(cfg,'chanstr'),         cfg.chanstr = 'noEyeABH';           end
if ~isfield(cfg,'slide'),           cfg.slide = 0.1;                    end
if ~isfield(cfg,'conds'),           cfg.conds = {'encLoc','encPer','encObj','encAni'};        end
if ~isfield(cfg,'subjplot'),        cfg.subjplot = [];                   end


doplots = cfg.doplots; adFile = cfg.adFile; chanstr = cfg.chanstr;
slide = cfg.slide; inconds = cfg.conds;

%load analysis details
load(adFile);

conds = {'encLoc','encPer','encObj','encAni'};
auc = nan(length(conds),1);
for icond = 1:length(conds)
    cond = conds{icond};
    if ismember(cond,inconds) %only open requested conds, but keep cond numbering consistent
        
        wname = sprintf('WRet_%s_%s_%0.2f.mat',conds{icond},chanstr,slide);
        wfname = fullfile(dirs.saveDirProc,substr,'ses1',wname);
        if ~exist(wfname,'file') error('W %s not found for subject %s\nlooking in %s',wname,substr,wfname);   end
        wfile.(cond) = load(wfname);
        
        trlord = ana.trl_order.(cond);
        if size(trlord,2)~=size(wfile.(cond).postrialinfo,2) %remove first column of trlord due to segmentation adding condnum
            trlord = cat(2,'condn',trlord);
            %warning('%s: adjusting trl_ord.%s to include condn column first!',substr,cond);
        end
        
        ordncol = ismember(trlord,'ordn');
        if sum(ordncol)~=1 error('too many sttp column found'); end
        crctcol = ismember(trlord,'crct');
        if sum(crctcol)~=1 error('nonunique crct column found'); end
        rtcol = ismember(trlord,'rtmm');
        if sum(rtcol)~=1 error('nonuniqe rt column found'); end
        
        pt = wfile.(cond).postrialinfo;
        nt = wfile.(cond).negtrialinfo;
        
        maxpos = max(wfile.(cond).avgtestsim,[],2);
        maxneg = max(wfile.(cond).negsim,[],2);
        p = cat(1,maxpos,maxneg);
        l = zeros(size(p)); l(1:length(maxpos)) = 1;
        [x,y,~,pvnauc] = perfcurve(l,p,1);
        
        if doplots
            if icond == 1 figure('color','white'); hold on; end
            plot(x,y,'linewidth',2);
        end
            
    end
    auc(icond) = pvnauc;    
end
if doplots
    legend(conds);
    title(substr);
end
wfile.auc = auc;

