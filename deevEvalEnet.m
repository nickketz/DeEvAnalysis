function out = deevEvalEnet(substr,cfg)
%function to select best performing elastic net per subject
% returns struct with predictions and evalutions 
% last dimension is 1:max diff 2:mean
%
% input;
%   substr: string; subject name to run
%   selmeas: string; eval measure to select model ('auc', def:'dprime')
%
% output:
%   out: struct with performance and selected predictions
%

if ~exist('cfg','var'),             cfg = [];                           end
if ~isfield(cfg,'selmeas')          cfg.selmeas = 'dprime';             end
if ~isfield(cfg,'slide')            cfg.slide = 2;                      end
if ~isfield(cfg,'tlen')             cfg.tlen = 5;                       end
if ~isfield(cfg,'doplots')          cfg.doplots = 0;                    end
if ~isfield(cfg,'adFile')           cfg.adFile = '/home/nike3851/work/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/v1/ft_data/OL_CL_encLoc_encPer_encObj_encAni_encNotLoc_encNotPer_encNotObj_encNotAni_eq0_art_ftAuto/pow_wavelet_w4_pow_3_50/analysisDetails.mat'; end


selmeas = cfg.selmeas; tslide = cfg.slide; tlen = cfg.tlen; doplots = cfg.doplots; adFile = cfg.adFile;

load(adFile);
savedir = fullfile(dirs.saveDirProc,substr,'ses1');
conds = {'encLoc','encPer','encObj','encAni'};


if doplots
    hfig = figure('color','white');
    c = get(gca,'colororder');
end
for icond = 1:length(conds)
    fprintf('\neval %s...',conds{icond});
    fname = sprintf('Enet_%s_%dtlen_%dtslide.mat',conds{icond},tlen,tslide);
    in = load(fullfile(savedir,fname));
        
    tmin = in.t(1); tmax = in.t(2); tdiff = .04; tdur = ((tmax-(tlen*tdiff))-tmin);
    nslide = size(in.fit,1);
    t = tmin:tslide*tdiff:tmax-tdiff;
    t = cat(1,t,t+tlen*tdiff);
    nconds = length(conds);
    if ~exist('auc','var')
        auc = nan(nslide,nconds,2); dp = nan(nslide,nconds,2); h = dp; fa = dp;%maxp x meanp
        b = cell(4,2); fit = b;
        tsel = nan(4,2,2); msel = nan(4,2);
    end

    
    l = in.test_lbl(:,icond);
    if icond == 1, p = nan([size(in.test_lbl),2]); end
    %maxp eval
    binvar = round(in.maxp);
    l = repmat(l,[1,size(binvar,2)]);
    h(:,icond,1) = sum(binvar==1 & l==1)./sum(l==1); 
    h(h==0) = .5/sum(l(:,1)==1); h(h==1) = 1-.5/sum(l(:,1)==1);
    fa(:,icond,1) = sum(binvar==1 & l==0)./sum(l==0); 
    fa(fa==0) = .5/sum(l(:,1)==0); fa(fa==1) = 1-.5/sum(l(:,1)==0);
    
    dp(:,icond,1) = norminv(h(:,icond,1)) - norminv(fa(:,icond,1));
    for islide = 1:nslide
        [~,~,~,auc(islide,icond,1)] = perfcurve(l(:,islide),in.maxp(:,islide),1);
    end
    %meanp eval
    binvar = round(in.meanp);
    h(:,icond,2) = sum(binvar==1 & l==1)./sum(l==1); h(h==0) = .5/sum(l(:,1)==1); h(h==1) = 1-.5/sum(l(:,1)==1);
    fa(:,icond,2) = sum(binvar==1 & l==0)./sum(l==0); fa(fa==0) = .5/sum(l(:,1)==0); fa(fa==1) = 1-.5/sum(l(:,1)==1);
    dp(:,icond,2) = norminv(h(:,icond,2)) - norminv(fa(:,icond,2));
    for islide = 1:nslide
        [~,~,~,auc(islide,icond,2)] = perfcurve(l(:,islide),in.meanp(:,islide),1);
    end
    fprintf('done');
    
    %select best for this condition
    for ieval = 1:2    
        switch selmeas
            case {'dprime'}
                [~,msel(icond,ieval)] = max(dp(:,icond,ieval));
            case {'auc'}
                [~,msel(icond,ieval)] = max(auc(:,icond,ieval));
            otherwise
                error('selection measure not recognized');
        end
                
        b{icond,ieval} = in.b{msel(icond,ieval)};
        fit{icond,ieval} = in.fit{msel(icond,ieval)};
        tsel(icond,:,ieval) = t(:,msel(icond,ieval));
        if ieval == 1, p(:,icond,ieval) = in.maxp(:,msel(icond,ieval)); end
        if ieval == 2, p(:,icond,ieval) = in.meanp(:,msel(icond,ieval)); end
    end

    if doplots
        subplot(2,2,icond);
        plot(mean(t),squeeze(auc(:,icond,:)),'-','linewidth',2);
        if icond == 1 legend('max','mean'); end            
        hold on
        set(gca,'colororderindex',1);
        plot(mean(t),squeeze(dp(:,icond,:)),'--','linewidth',2);
        title(conds{icond});
        box off
    end
            
end
eval.auc = auc;
eval.dprime = dp;
eval.t = t;
eval.hit = h;
eval.fa = fa;

out.pred = p;
out.label = in.test_lbl;
out.b = b;
out.fit = fit;
out.tsel = tsel;
out.msel = msel;
out.t = t;
out.trialinfo = in.test_trialinfo;
out.eval = eval;
out.selmeas = selmeas;

% figure('color','white')
% subplot(2,2,1);
% plot(mean(t),auc(:,:,1),'linewidth',2);
% legend(conds);
% title('maxp auc');
% box off
% subplot(2,2,3);
% plot(mean(t),dp(:,:,1),'linewidth',2);
% title('maxp dprime');
% box off
% subplot(2,2,2);
% plot(mean(t),auc(:,:,2),'linewidth',2);
% title('meanp auc');
% box off
% subplot(2,2,4);
% plot(mean(t),dp(:,:,2),'linewidth',2);
% title('meanp dprime');
% box off

    
    