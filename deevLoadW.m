function [wfile,subjplot] = deevLoadW(substr,cfg)
%
% function to get weight vector W for a given subject across encoding
% conditions, also grabs performance of W across OL and CL, and correct vs.
% incorrct retrievals
%
% input:
%   substr: string of subject to retrieve W info for
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

if doplots
    if isempty(cfg.subjplot)
        hdiff = figure('color','white','Name','OL vs. CL');
        hcrct = figure('color','white','name','Corrct vs. Incorrect');
    else
        hdiff = cfg.subjplot(1);
        hcrct = cfg.subjplot(2);
    end
end

conds = {'encLoc','encPer','encObj','encAni'};

for icond = 1:length(conds)
    cond = conds{icond};
    if ismember(cond,inconds) %only open requested conds, but keep cond numbering consistent
        
        wname = sprintf('W_%s_%s_%0.2f.mat',conds{icond},chanstr,slide);
        wfname = fullfile(dirs.saveDirProc,substr,'ses1',wname);
        if ~exist(wfname,'file') error('W %s not found for subject %s\nlooking in %s',wname,substr,wfname);   end
        wfile.(cond) = load(wfname);
        
        trlord = ana.trl_order.(cond);
        if size(trlord,2)~=size(wfile.(cond).cltrialinfo,2) %remove first column of trlord due to segmentation adding condnum
            trlord = cat(2,'condn',trlord);
            warning('%s: adjusting trl_ord.%s to include condn column first!',substr,cond);
        end
        
        ordncol = ismember(trlord,'ordn');
        if sum(ordncol)~=1 error('too many sttp column found'); end
        crctcol = ismember(trlord,'crct');
        if sum(crctcol)~=1 error('nonunique crct column found'); end
        rtcol = ismember(trlord,'rtmm');
        if sum(rtcol)~=1 error('nonunqie rt column found'); end
        
        %find correct trials
        clcrct = wfile.(cond).cltrialinfo(:,crctcol);
        olcrct = wfile.(cond).oltrialinfo(:,crctcol);
        
        %find trials where this element is the target
        clsttp = cellfun(@num2str,num2cell(wfile.(cond).cltrialinfo(:,ordncol)),'uniformoutput',0);
        clretind = ~cellfun(@isempty,regexp(clsttp,[num2str(icond) '$'])); %trials that include the current condition
        olsttp = cellfun(@num2str,num2cell(wfile.(cond).oltrialinfo(:,ordncol)),'uniformoutput',0);
        olretind = ~cellfun(@isempty,regexp(olsttp,[num2str(icond) '$']));
        
        %find trials where this element is the non-target
        
        %compare ol and cl on retrieval trials of this condition
        wfile.(cond).olscore_crct = max(wfile.(cond).olsim(olretind & olcrct,:),[],2);
        wfile.(cond).olscore_incrct = max(wfile.(cond).olsim(olretind & ~olcrct,:),[],2);
        
        wfile.(cond).clscore_crct = max(wfile.(cond).clsim(clretind & clcrct,:),[],2);
        wfile.(cond).clscore_incrct = max(wfile.(cond).clsim(clretind & ~clcrct,:),[],2);
        
        wfile.(cond).olscore = max(wfile.(cond).olsim(olretind,:),[],2);
        wfile.(cond).clscore = max(wfile.(cond).clsim(clretind,:),[],2);
        
        if doplots
            e = 0:.01:1;
            figure(hdiff);
            a = subplot(2,2,icond);
            hold on
            corder = get(gca,'colororder');
            histogram(a,wfile.(cond).clscore, e, 'normalization','probability');
            histogram(a,wfile.(cond).olscore, e, 'normalization','probability');
            legend({'CL','OL'});
            plot(ones(2,1)*mean(wfile.(cond).clscore), ylim, '--','color',corder(1,:),'linewidth',3);
            plot(ones(2,1)*mean(wfile.(cond).olscore), ylim, '--','color',corder(2,:),'linewidth',3);
            title([substr '-' cond]);
            
            e = 0:.01:1;
            figure(hcrct);
            a = subplot(2,2,icond);
            hold on
            crct = cat(1,wfile.(cond).clscore_crct,wfile.(cond).olscore_crct);
            incrct = cat(1,wfile.(cond).clscore_incrct,wfile.(cond).olscore_incrct);
            corder = get(gca,'colororder');
            histogram(a,crct, e, 'normalization','probability');
            histogram(a,incrct, e, 'normalization','probability');
            legend({'crct','incrct'});
            plot(ones(2,1)*mean(crct), ylim, '--','color',corder(1,:),'linewidth',4);
            plot(ones(2,1)*mean(incrct), ylim, '--','color',corder(2,:),'linewidth',4);
            title([substr '-' cond]);
        end
    end
    
end
subjplot = [];
if cfg.doplots   
    subjplot(1) = hdiff;
    subjplot(2) = hcrct;    
end

