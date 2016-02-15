function deev_seg_wrapper(whichStages)
% tnt_ftprocess_tfr_pow_wrapper(whichStages)
%
% To run on dream, at the command line type: distmsub core_conn_wrapper.m
%
% To run on a local computer, type the command in MATLAB
%
% There is only one stage:
%  stage1 = call wrapper that calls create_ft_struct (which calls seg2ft,
%  which calls ft_freqanalysis) and saves one file per subject
%
% Input:
%  whichStages: the stage number(s) to run (default = 1)
%
% Output:
%  time-frequency data

% add eeg_directories to the path
%add_eeg_path

% check/handle arguments
error(nargchk(0,1,nargin))
STAGES = 1;
if nargin == 1
    STAGES = whichStages;
end

runLocally = 1;

% load analysisDetails instead of using specifications below
%if runLocally == 0
adFile = '/projects/nike3851/CORE_CONN/eeg/1000to2000/ft_data/Word_Blue_Face_Word_Blue_Scene_Word_Green_Face_Word_Green_Scene_Word_Red_Face_Word_Red_Scene_eq0_art_nsAuto/conn_wavelet_w4_fourier_-500_1980_3_50/analysisDetails.mat';

[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile);

%  elseif runLocally == 1
%    adFile = '/Volumes/curranlab/TNT/TNT_matt/eeg/-1000_1700/ft_data/B_NT_TH_eq1/conn_scd_mtmconvol_hanning_fourier_-500_980_3_9/analysisDetails.mat';
%    [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,1);
%end

%% set up for running stages and specifics for Dream

% name(s) of the functions for different stages of processing
stageFun = {@stage1};
timeOut  = {12}; % in HOURS

if runLocally == 0
    % need to export DISPLAY to an offscreen buffer for MATLAB DCS graphics
    sched = findResource();
    if strcmp(sched.Type, 'generic')
        setenv('DISPLAY', 'dream:99');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%capture diary and time statistics
thisRun = [exper.name,'_overview_',datestr(now,'ddmmmyyyy-HHMMSS')];
%thisRun = [exper.name,'_overview_',datestr(now,7) datestr(now,3) datestr(now,10)];
diary(fullfile(dirs.saveDirProc,[thisRun '.log']));
tStart = tic;
fprintf('START TIME: %s\n',datestr(now,13));
for i = STAGES
    tS = tic;
    fprintf('STAGE%d START TIME: %s\n',i, datestr(now,13));
    
    % execute the processing stage
    stageFun{i}(ana,exper,dirs,runLocally,timeOut{i});
    
    fprintf('STAGE%d END TIME: %s\n',i, datestr(now,13));
    fprintf('%.3f -- elapsed time STAGE%d (seconds)\n', toc(tS), i);
end
time = toc(tStart);
fprintf('%.3f -- elapsed time OVERALL (seconds)\n', time);
fprintf('END TIME: %s\n',datestr(now,13));
diary off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stage1(ana,exper,dirs,runLocally,timeOut)
% stage1: process the input files with FieldTrip based on the analysis
% parameters

%% Process the data
if runLocally == 0
    %% Dream: create one task for each subject (i.e., submit one to each node)
    
    % start a new job
    job = newJob(dirs);
    
    cfg_ft = [];
    cfg_ft.method = 'wpli_debiased';
    %  cfg_ft.method = 'coh';
    %cfg_ft.method = 'plv';
    
    clusnames = fieldnames(ana.stat_clus);
    for sub = 1:1%length(exper.subjects)
        for ses = 1:length(exper.sessions)
            for clus = 1:length(clusnames)
                
                %get cluster electrodes as cell array of elec names
                mask = ana.stat_clus.(clusnames{clus}).mask;
                mask = sum(sum(mask==1,3),2)>0;
                elecs = ana.stat_clus.(clusnames{clus}).label;
                cluselec = elecs(mask);
                notclus = elecs(~mask);
                
                % do all pairwise combinations with elec outside of cluster
                chans_to_pair = cell(length(notclus)*length(cluselec),2);
                ind = 0;
                for i = 1:length(cluselec)
                    for j = 1:length(notclus)
                        ind = ind+1;
                        chans_to_pair{ind,1} = cluselec{i};
                        chans_to_pair{ind,2} = notclus{j};
                    end
                end
                cfg_ft.channelcmb = chans_to_pair;
                for evVal = 2:length(ana.eventValues{1})
                    
                    cfg_ft.inputfile = fullfile(dirs.saveDirProc,exper.subjects{sub},sprintf('ses%d',ses),sprintf('data_%s_%s.mat',ana.ftype,ana.eventValues{1}{evVal}));
                    cfg_ft.outputfile = fullfile(dirs.saveDirProc,exper.subjects{sub},sprintf('ses%d',ses),sprintf('data_%s_%s_%s.mat',cfg_ft.method,ana.eventValues{1}{evVal},clusnames{clus}));
                    
                    
                    %if ~exist(cfg_ft.outputfile,'file')
                    %  fprintf('Processing %s, %s, %s...\n',exper.subjects{sub},exper.sessions{ses},exper.eventValues{evVal});
                    ft_connectivityanalysis(cfg_ft);
                    fprintf('Done.\n');
                    %else
                    %  fprintf('ALREADY EXISTS: %s\n',cfg_ft.outputfile);
                    %end
                end
            end
        end
    end
    
    
    runJob(job, timeOut, fullfile(dirs.saveDirProc,[exper.name,'_stage1_',datestr(now,'ddmmmyyyy-HHMMSS'),'.log']));
    
    % final step: destroy the job because this doesn't happen in runJob
    destroy(job);
    
else
    %% run the function locally
    
    % create a log of the command window output
    thisRun = [exper.name,'_stage1_',datestr(now,'ddmmmyyyy-HHMMSS')];
    % turn the diary on
    diary(fullfile(dirs.saveDirProc,[thisRun,'.log']));
    
    % use the peer toolbox
    %ana.usePeer = 1;
    ana.usePeer = 0;
    
    cfg_ft = [];
    cfg_ft.method = 'wpli';
    
    %   % do some pairwise combinations
    %   chans_to_pair = {'E1'; 'E6'; 'E9'; 'E11'; 'E15'; 'E22'; 'E24'; 'E29'; 'E32'; 'E33'; 'E34'; 'E37'; 'E41'; 'E43'; 'E45'; 'E51'; 'E52'; 'E53'; 'E55'; 'E62'; 'E64'; 'E66'; 'E69'; 'E72'; 'E75'; 'E81'; 'E84'; 'E86'; 'E87'; 'E89'; 'E92'; 'E95'; 'E97'; 'E103'; 'E108'; 'E111'; 'E116'; 'E120'; 'E122'; 'E124'};
    %   cfg_ft.channelcmb = nchoosek(chans_to_pair,2);
    
    
    clusnames = fieldnames(ana.stat_clus);
    for sub = 1:1%length(exper.subjects)
        for ses = 1:length(exper.sessions)
            for clus = 1:length(clusnames)
                
                %get cluster electrodes as cell array of elec names                
                mask = ana.stat_clus.(clusnames{clus}).mask;
                mask = sum(sum(mask==1,3),2)>0;
                elecs = ana.stat_clus.(clusnames{clus}).label;
                cluselec = elecs(mask);
                notclus = elecs(~mask);
                
                % do all pairwise combinations with elec outside of cluster
                chans_to_pair = cell(length(notclus)*length(cluselec),2);
                ind = 0;
                for i = 1:length(cluselec)
                    for j = 1:length(notclus)
                        ind = ind+1;
                        chans_to_pair{ind,1} = cluselec{i};
                        chans_to_pair{ind,2} = notclus{j};
                    end
                end
                cfg_ft.channelcmb = chans_to_pair;
                %cfg_ft.channelcmb = {'all','all'};
                for evVal = 1:length(ana.eventValues{1})
                    
                    cfg_ft.inputfile = fullfile(dirs.saveDirProc,exper.subjects{sub},sprintf('ses%d',ses),sprintf('data_%s_%s.mat',ana.ftype,ana.eventValues{1}{evVal}));
                    cfg_ft.outputfile = fullfile(dirs.saveDirProc,exper.subjects{sub},sprintf('ses%d',ses),sprintf('data_%s_%s_%s.mat',cfg_ft.method,ana.eventValues{1}{evVal},clusnames{clus}));
                    
                    
                    %if ~exist(cfg_ft.outputfile,'file')
                    %  fprintf('Processing %s, %s, %s...\n',exper.subjects{sub},exper.sessions{ses},exper.eventValues{evVal});
                    ft_connectivityanalysis(cfg_ft);
                    fprintf('Done.\n');
                    %else
                    %  fprintf('ALREADY EXISTS: %s\n',cfg_ft.outputfile);
                    %end
                end
            end
        end
    end
    
    % turn the diary off
    diary off
end

function job = newJob(dirs)
% newJob Creates a new PCT job and sets job's dependencies
%
%   dirs -- data structure with necessary fields like data locations

% Set up scheduler, job
sched = findResource();
job = createJob(sched);
% define the directories to add to worker sessions' matlab path
homeDir = getenv('HOME');
myMatlabDir = fullfile(homeDir,'Documents','MATLAB');
p = path();
set(job, 'PathDependencies', {homeDir, myMatlabDir, pwd(), p, dirs.dataroot});

function runJob( job, timeOut, logFile )
% runJob Submits and waits on job to finish or timeout
%  runJob will submit the supplied job to the scheduler and will
% wait for the job to finish or until the timeout has been reached.
% If the job finishes, then the command window outputs of all tasks
% are appended to the log file and the job is destroyed.
%   If the timeout is reached, an error is reported but the job is not
% destroyed.
%
%   job -- the job object to submit
%   timeOut -- the timeout value in hours
%   logFile -- full file name of the log file to append output to
%
% Example:
%       runJob( job, 5, 'thisrun.log');

% check/handle arguments
error(nargchk(1,3,nargin))
TIMEOUT=3600*5; % default to 5 hours
if nargin > 1
    TIMEOUT=timeOut*3600;
end
LOGFILE=[job.Name '.log'];
if nargin > 2
    LOGFILE = logFile;
end

% Capture command window output from all tasks
alltasks = get(job, 'Tasks');
set(alltasks, 'CaptureCommandWindowOutput', true);

% Submit Job/Tasks and wait for completion (or timeout)
submit(job)
finished = waitForState(job, 'finished', TIMEOUT);
if finished
    errors = logOutput(alltasks, LOGFILE);
    if errors
        error([mfilename ':logOutput'],'%s had %d errors',job.Name, errors)
        %elseif ~errors
        %  destroy(job);
    end
else
    error([mfilename ':JobTimeout'],'%s: Timed out waiting for job...NAME: %s',...
        datestr(now, 13), job.Name, job.ID, job.StartTime)
end

function numErrors=logOutput( tasks, logFile )
% logOutput - concatenates tasks' output into a logfile
%   tasks -- the tasks to capture output from
%   logFile -- the file to log the output to
%   numErrors -- number of tasks which failed

% check for argument(s)
error(nargchk(2,2,nargin))

numErrors=0;
try
    fid=fopen(logFile, 'a+');
    for i=1:length(tasks)
        fprintf(fid,'\n***** START TASK %d *****\n',i);
        fprintf(fid,'%s\n', tasks(i).CommandWindowOutput);
        if ~isempty(tasks(i).Error.stack)
            numErrors = numErrors +1;
            % write to log file
            fprintf( fid, 'ERROR: %s\n', tasks(i).Error.message );
            fprintf( fid, '%s\n', tasks(i).Error.getReport );
            % write to standard error
            fprintf( 2, 'ERROR: %s\n', tasks(i).Error.message );
            fprintf( 2, '%s\n', tasks(i).Error.getReport );
        end
        fprintf(fid,'\n***** END TASK %d *****\n',i);
    end
    fclose(fid);
catch ME
    disp(ME)
    warning([mfilename ':FailOpenLogFile'],...
        'Unable to write log file with task output...');
end

% function add_eeg_path
%
% %% add specific toolboxes to the path - you could instead put this in ~/Documents/MATLAB/startup.m
%
% %% initialize
% homeDir = getenv('HOME');
%
% myMatlabDir = fullfile(homeDir,'Documents','MATLAB');
%
% %% set up eegToolbox path
% eegToolboxDir = dir(fullfile(myMatlabDir,'eegToolbox*'));
% if ~isempty(eegToolboxDir)
%   eegToolboxDir = fullfile(myMatlabDir,eegToolboxDir.name);
%   % add top folder and all subfolders
%   addpath(genpath(eegToolboxDir));
% end
%
% %% set up eeglab path
% eeglabDir = dir(fullfile(myMatlabDir,'eeglab*'));
% if ~isempty(eeglabDir)
%   eeglabDir = fullfile(myMatlabDir,eeglabDir.name);
%   % add top folder and all subfolders
%   addpath(genpath(eeglabDir));
%
%   % remove eeglab's external directory if it was added
%   eeglabExtDir = fullfile(eeglabDir,'external');
%   if ~isempty(eeglabExtDir)
%     rmpath(genpath(eeglabExtDir));
%   end
%
%   % % remove eeglab's fieldtrip directory if it was added
%   % eeglabFtDir = dir(fullfile(eeglabDir,'external','fieldtrip*'));
%   % if ~isempty(eeglabFtDir)
%   %   eeglabFtDir = fullfile(myMatlabDir,eeglabFtDir.name);
%   %   rmpath(genpath(eeglabFtDir));
%   % end
% end
%
% %% set up fieldtrip path
% ftDir = dir(fullfile(myMatlabDir,'fieldtrip*'));
% if ~isempty(ftDir)
%   ftDir = fullfile(myMatlabDir,ftDir.name);
%   % add only the top folder
%   addpath(ftDir);
%   % run ft_defaults to add the subdirectories that FT needs
%   ft_defaults
%
%   % add the peer directory
%   addpath(fullfile(ftDir,'peer'));
%
%   % add the SPM directory
%   %addpath(fullfile(ftDir,'external','spm8'));
%
%   % % remove fieldtrip's external directory
%   % ftExtDir = fullfile(ftDir,'external');
%   % if ~isempty(ftExtDir)
%   %   rmpath(genpath(ftExtDir));
%   % end
%
%   % remove fieldtrip's eeglab directory if it was added
%   ftEeglabDir = dir(fullfile(ftDir,'external','eeglab*'));
%   if ~isempty(ftEeglabDir)
%     ftEeglabDir = fullfile(myMatlabDir,ftEeglabDir.name);
%     rmpath(genpath(ftEeglabDir));
%   end
% end
%
% %% set up EP_Toolkit path
% epDir = dir(fullfile(myMatlabDir,'EP_Toolkit*'));
% if ~isempty(epDir)
%   epDir = fullfile(myMatlabDir,epDir.name);
%   % add top folder and all subfolders
%   addpath(genpath(epDir));
% end
%
% %% add my other analysis scripts/files
% addpath(genpath(fullfile(myMatlabDir,'recogmodel_mvm')));
%
% %% add my experiment, fieldtrip, and RM ANOVA scripts
% addpath(genpath(fullfile(myMatlabDir,'mat_mvm')));
%
% %% put the ~/Documents/MATLAB folder at the top of the path
% addpath(myMatlabDir);
%
% %% remove CVS and .svn directories from path
% entries = regexp(path, ['[^',pathsep,']*',pathsep], 'match');
% for i = 1:length(entries)
%   entry = char(entries{i});
%   if ~isempty(strfind(entry, '.svn'))
%     rmpath(entry);
%   end
%   if ~isempty(strfind(entry, 'CVS'))
%     rmpath(entry);
%   end
% end
%
