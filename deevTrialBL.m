function [data,exper] = deevTrialBL(exper,dirs,ana,cfg)
% trial based baseline correction within subject
% do baseline correction of data.(cfg.trgCond) with matching trials from cfg.twin of
% data.(cfg.blCond).  blCond and trgCond musth have same trialinfo trl_order
%
% input:
%   exper: experiment structure
%   dirs: dirs structure
%   ana: ana structure configed to load blCond and trgCond, this can use
%        trialinfo loading
%   cfg: config struct
%       .param = str; name of field with the data (e.g. 'powspctrm')
%       .loadcfg = cfg struct to pass into mm_ft_loadData_multises
%       .colind = 1xN logical; N is size(ana.trl_order.(trgCond),2),
%          columns that must match for trials to be considered the sams
%       .twin = 1x2 float; time window to create average basline 
%       .blCond = str; condition to use as basline
%       .trgCond = str; target condition to get corrected
%       .method = str; zscore(def),absolute,relative,relchange
%       .keeptrials = str; 'yes' outputs trials after bl correction, 'no' averages them
%
% output:
%   data: data struct with new condition bl(cfg.trgcond)
%

if ~isfield(cfg,'method') cfg.method = 'zscore';        end

%must keeptrials
cfg.loadcfg.keeptrials = 'yes';
data = [];

fprintf('\n%s baseline correcting %s in reference to %s...',cfg.method,cfg.trgCond,cfg.blCond);
for isub = 1:length(exper.subjects)
    
    try
        %load trgcond and blcond
        subexper = nk_rmSubs(exper,~strcmp(exper.subjects{isub},exper.subjects));
        indata = mm_ft_loadData_multiSes(cfg.loadcfg,subexper,dirs,ana);
    catch
        %load data failed return; empty slot in data for this sub
        data.sub(isub).data = [];
        continue
    end
    
    if isfield(ana,'eventValuesComb')
        cfg_appnd = []; cfg_appnd.parameter = cfg.param; cfg_appnd.appenddim = 'rpt';
        for ises = 1:length(ana.eventValuesComb)
            for inew = 1:length(ana.eventValuesNew{ises})
                cmbstr = '';
                for icmb = 1:length(ana.eventValuesComb{ises}{inew})
                    cmbstr = sprintf('%s, indata.%s.%s.sub.data',cmbstr,exper.sesStr{1},ana.eventValuesComb{ises}{inew}{icmb});
                end
                eval(['cmbdata = ft_appendfreq(cfg_appnd' cmbstr ');']);
                indata.(exper.sesStr{1}) = rmfield(indata.(exper.sesStr{1}),ana.eventValuesComb{ises}{inew});
                indata.(exper.sesStr{1}).(ana.eventValuesNew{ises}{inew}).sub.data = cmbdata;
            end
        end
    end
                

    bldata = indata.(exper.sesStr{1}).(cfg.blCond).sub.data;
    trgdata = indata.(exper.sesStr{1}).(cfg.trgCond).sub.data;
    clear indata cmbdata;
    
    %find matching trials and sort bl to correspond to trg
    [blintrg, trginbl] = ismember(bldata.trialinfo(:,cfg.colinds),trgdata.trialinfo(:,cfg.colinds),'rows');
    trginbl = trginbl(blintrg);
    %subselect bldata
    bldata.(cfg.param) = bldata.(cfg.param)(blintrg,:,:,:); 
    bldata.trialinfo = bldata.trialinfo(blintrg,:);
    %subselect trgdata and subfields
    trgdata.(cfg.param) = trgdata.(cfg.param)(trginbl,:,:,:);
    trgdata.trialinfo = trgdata.trialinfo(trginbl,:);
    if isfield(bldata,'cumtapcnt')
        trgdata.cumtapcnt = trgdata.cumtapcnt(trginbl,:);
    end
    %sort bl
    [blintrg, trginbl] = ismember(bldata.trialinfo(:,cfg.colinds),trgdata.trialinfo(:,cfg.colinds),'rows');
    if any(blintrg~=1) || any(trginbl==0), error('sub selection has gone wrong'); end
    bldata.(cfg.param) = bldata.(cfg.param)(trginbl,:,:,:);
    bldata.trialinfo = bldata.trialinfo(trginbl,:);
    if isfield(bldata,'cumtapcnt'), bldata.cumtapcnt = bldata.cumtapcnt(trginbl,:); end
    
    if sum(sum(bldata.trialinfo(:,cfg.colinds) - trgdata.trialinfo(:,cfg.colinds))) > 0
        error('sorting baseline data has gone wrong');
    end
    
    %set sub data struct
    data.sub(isub).data = trgdata;
    
    %get bl average
    blt = bldata.time>=cfg.twin(1) & bldata.time<=cfg.twin(2);
    if ndims(bldata.(cfg.param))~=4 error('only works for non-trial averaged data'); end
    blm = nanmean(bldata.(cfg.param)(:,:,:,blt),4);
    
    if strcmp(cfg.method,'zscore')
        fprintf('Zscoring data relative to %s mean([%.2f %.2f].\n',cfg.blCond,cfg.twin(1),cfg.twin(2));
        % std across time, then avg across events (lower freqs often get smaller std)
        blstd = nanmean(nanstd(trgdata.(cfg.param)(:,:,:,blt),0,4),1);
        trgdata.(cfg.param) = bsxfun(@rdivide,bsxfun(@minus,trgdata.(cfg.param),blm),blstd);
    elseif strcmp(cfg.method,'absolute')
        fprintf('Subtracting %s mean([%.2f %.2f]) power from entire trial.\n',cfg.blCond,cfg.twin(1),cfg.twin(2));
        trgdata.(cfg.param) = bsxfun(@minus,trgdata.(cfg.param),blm);
    elseif strcmp(cfg.method,'relative')
        fprintf('Dividing entire trial power by %s mean([%.2f %.2f]) power.\n',cfg.blCond,cfg.twin(1),cfg.twin(2));
        trgdata.(cfg.param) = bsxfun(@rdivide,trgdata.(cfg.param),blm);
    elseif strcmp(cfg.method,'relchange')
        fprintf('Subtracting %s mean([%.2f %.2f]) power from entire trial and dividing by that mean.\n',cfg.blCond,cfg.twin(1),cfg.twin(2));
        trgdata.(cfg.param) = bsxfun(@rdivide,bsxfun(@minus,trgdata.(cfg.param),blm),blm);
    end
    
    
    %null data after response
    if cfg.nullRT
        fprintf('\nnulling data after reponse time...');
        
        if ~isfield(cfg,'bufRT'), cfg.bufRT = 0; end %subtract buffer from RT end? i.e. for button press/stim onset
        %get rt data
        rtcol = ismember(ana.trl_order.(ana.eventValues{1}{1}),'rtmm');
        if sum(rtcol)~=1, error('non-unique column found for rtmm'); end
        rt = trgdata.trialinfo(:,rtcol) - cfg.bufRT;
        %is this response locked or stimulus locked?
        isRL = 0;
        if sum(trgdata.time<0) > sum(trgdata.time>0), rt = rt*-1; isRL=1;  end

        for itrl = 1:size(trgdata.trialinfo,1)
            %find time points to nan for this trial
            if isRL
                tsel = trgdata.time<=rt(itrl);
            else
                tsel = trgdata.time>=rt(itrl);
            end
            if sum(tsel)>0, trgdata.(cfg.param)(itrl,:,:,tsel) = nan;   end
        end
    end
    
    %avg over trials?
    if strcmp('no',cfg.keeptrials); 
        fprintf('averaging over trials\n');
        trgdata = ft_selectdata(struct('avgoverrpt','yes','nanmean','yes'),trgdata);
    end
    
    % put in the trial counts
    exper.nTrials.(['tbl' cfg.trgCond])(isub,1) = size(bldata.(cfg.param),1);
    fprintf('\n---\n%s %s returned %d trials\n---\n\n',exper.subjects{isub},cfg.trgCond,size(bldata.(cfg.param),1));
    
    %append to data struct
    data.sub(isub).data = trgdata;
    
end

fprintf('done\n');


