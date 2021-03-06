function [comp,cfg_sasica] = deevSASICA(comp,data,ana,dirs)

addpath(fullfile(dirs.dataroot,'..','scripts','SASICA'));
cfg_sasica = [];

cfg_sasica.opts.plotallcomp = 0;

cfg_sasica.chancorr.enable = 0;
cfg_sasica.autocorr.enable = 1;
cfg_sasica.EOGcorr.enable = 1;
cfg_sasica.SNR.enable = 0;
cfg_sasica.trialfoc.enable = 0;

cfg_sasica.FASTER.enable = 1;
cfg_sasica.ADJUST.enable = 0;
cfg_sasica.MARA.enable = 0;

if isfield(comp,'reject')
    %sasica already run
    cfg_sasica.opts.nocompute = 1;
    fprintf('not redoing computations, just replotting\n');
else
    
    ana = mm_ft_elecGroups(ana);
    heye = ana.elecGroups{strcmp('eyes_horiz',ana.elecGroupsStr)};
    aeye = ana.elecGroups{strcmp('eyes_above',ana.elecGroupsStr)};
    beye = ana.elecGroups{strcmp('eyes_below',ana.elecGroupsStr)};
    veye = cat(2,aeye,beye);
    
    ts = cat(2,data.trial{:});
    VEOG = ts(ismember(data.label,aeye),:) - ts(ismember(data.label,beye),:);
    cfg_sasica.EOGcorr.VEOG = mean(VEOG,1);
    cfg_sasica.EOGcorr.HEOG = ts(ismember(data.label,heye{1}),:) - ts(ismember(data.label,heye{2}),:);
        
    cfg_sasica.FASTER.blinkchans = cat(2,veye,heye);
        
    comp.chanlocs = readlocs('GSN-HydroCel-129.sfp');
    comp.chanlocs = comp.chanlocs(4:end-1); %remove fid chanels and CZ
    comp.layout = ft_prepare_layout(struct('elecfile','GSN-HydroCel-129.sfp'));
    comp.data = data.trial;
    %comp.icasphere = 2.0*inv(sqrtm(cov(ts(1:end-1,:)')));
end

try
    evalin('base','clear comp');  %make sure there's no preexisting comp in base workspace
    ft_SASICA(comp,cfg_sasica);
    uiwait(gcf);
    %grab comp from base, then remove it 
    comp = evalin('base','comp');
    evalin('base','clear comp');     
catch
    warning('error in plotting, redoing...');
end
rmpath(fullfile(dirs.dataroot,'..','scripts','SASICA'));

