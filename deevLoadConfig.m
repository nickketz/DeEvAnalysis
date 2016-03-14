function cfg = deevLoadConfig()
% function to return default config struct for loading pow data
%
% input:
%
% output: 
%   cfg: config struct with default parameters for loading subject power
%   data
%


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
