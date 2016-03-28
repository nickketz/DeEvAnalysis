function data = deevCondBL(data,cfg)
% condition based baseline correction within subject
% do baseline correction of data.(cfg.trgCond) with time average from cfg.twin of
% data.(cfg.blCond). Adds within subject baseline corrected condition bl(cfg.trgCond) to data_pow
%
% input:
%   data: mm_ft freq data struct
%   cfg: config struct
%       .param = str; name of field with the data (e.g. 'powspctrm')
%       .twin = 1x2 vec; time window to average over
%       .blCond = str; condition to use as basline
%       .trgCond = str; target condition to get corrected
%       .method = str; zscore(def),absolute,relative,relchange
%
% output:
%   data: data struct with new condition bl(cfg.trgcond)
%

if ~isfield(cfg,'method') cfg.method = 'zscore';        end

ses = fieldnames(data);
for ises = 1:length(ses)
    %copy to new baseline corrected condition
    tmpcond = data.(ses{ises}).(cfg.trgCond);
    fprintf('\n%s baseline correcting %s in reference to %s...',cfg.method,cfg.trgCond,cfg.blCond);
    for isub = 1:length(tmpcond.sub)
        %get bl average
        bldata = data.(ses{ises}).(cfg.blCond).sub(isub).data;        
        blt = bldata.time>=cfg.twin(1) & bldata.time<=cfg.twin(2); 
        if ndims(bldata.(cfg.param))~=3 error('only works for trial averaged data'); end
        blm = nanmean(bldata.(cfg.param)(:,:,blt),3);
        
        if strcmp(cfg.method,'zscore')
            blstd = nanstd(bldata.(cfg.param)(:,:,blt),0,3);         
            tmpcond.sub(isub).data.(cfg.param) = bsxfun(@rdivide,bsxfun(@minus,tmpcond.sub(isub).data.(cfg.param),blm),blstd);
        elseif strcmp(cfg.method,'absolute')
            %fprintf('\tSubtracting mean([%.2f %.2f]) power from entire trial.\n',cfg.baseline_time(1),cfg.baseline_time(2));
            tmpcond.sub(isub).data.(cfg.param) = bsxfun(@minus,tmpcond.sub(isub).data.(cfg.param),blm);            
        elseif strcmp(cfg.method,'relative')
            %fprintf('\tDividing entire trial power by mean([%.2f %.2f]) power.\n',cfg.baseline_time(1),cfg.baseline_time(2));
            tmpcond.sub(isub).data.(cfg.param) = bsxfun(@rdivide,tmpcond.sub(isub).data.(cfg.param),blm);            
        elseif strcmp(cfg.method,'relchange')
            %fprintf('\tSubtracting mean([%.2f %.2f]) power from entire trial and dividing by that mean.\n',cfg.baseline_time(1),cfg.baseline_time(2));
            tmpcond.sub(isub).data.(cfg.param) = bsxfun(@rdivide,bsxfun(@minus,tmpcond.sub(isub).data.(cfg.param),blm),blm);
        end
    end 
    eval(sprintf('data.(ses{ises}).bl%s = tmpcond;',cfg.trgCond));
    fprintf('done\n');
end

