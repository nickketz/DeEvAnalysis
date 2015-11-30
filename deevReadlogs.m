function outdata = deevReadlogs(cfg)
% DeEv analysis script that reads in deev_*.txt files in 'logs/'
%
%   input: 
%       cfg
%        filefiltstr: regexp filter on file listing, default is 'EntAssoc_*_Sub[0-9]{1,2}\.txt'
%        badsubs: subject numbers to exclude from analysis
%        dir: struct containing the logs (def: logfiles)
%
%   output:
%       outdata: struc with headernames as matricies with Subs on the
%       columns




%defaults
if ~exist('cfg','var')          cfg = [];                end
if ~isfield(cfg,'badsubs')      cfg.badsubs = [];        end %which subjects to remove
if ~isfield(cfg,'filefiltstr')  cfg.filefiltstr = 'deevBehExp_36events_Sub[0-9]{1,2}\.txt';     end
if ~isfield(cfg,'dir')          cfg.dir = 'logfiles';    end %where are the log files?


files = dir([cfg.dir filesep '*.txt']);
files = {files.name};

filesmatch = regexp(files,cfg.filefiltstr);
filesmatch = ~cellfun(@isempty,filesmatch);
if sum(filesmatch)==0
    error('no log files found');
end
files = files(filesmatch);
for ifile = 1:length(files)
    files{ifile} = [cfg.dir filesep files{ifile}];
end

%remove bad subs
snums = regexp(files,'_Sub([0-9]+)\.txt','tokens');
snums = cellfun(@(x) (str2num(x{1}{1})),snums);
files = files(~ismember(snums,cfg.badsubs));
goodsubs = regexp(files,'_Sub([0-9]+)\.txt','tokens');
goodsubs = cellfun(@(x) (str2num(x{1}{1})),goodsubs);


%get var names and types from header
temp = textread(files{1},'%s');
vars = strread(temp{1},'%s','delimiter',';');
fstring = [];
vartype = {};
for ivar = 1:length(vars)
    switch vars{ivar}(end)
        case '$'
            fstring = [fstring '%s '];
            %outdata.(vars{ivar}(1:end-1)) = {};
            vartype{ivar} = '$';
        case '|'
            fstring = [fstring '%d '];
            %outdata.(vars{ivar}(1:end-1)) = [];
            vartype{ivar} = '|';
        case '#'
            fstring = [fstring '%f '];
            %outdata.(vars{ivar}(1:end-1)) = [];
            vartype{ivar} = '#';
        otherwise 
            fstring = [fstring '%s '];
            %outdata.(vars{ivar}(1:end-1)) = {};
            vartype{ivar} = '$';
    end
    vars{ivar} = vars{ivar}(1:end-1);    
end

%get data from logs
%mydata = cell(1,filesmatch);
cdata = cell(1,length(vars)+1);
for ilog = 1:length(files)    
    %read in log file
    logfilename = files{ilog};
    fid = fopen(logfilename,'r');
    mydata = textscan(fid, fstring, 'Headerlines', 1, 'Delimiter', ';', 'TreatAsEmpty', {'na'});        
    fclose(fid);
    
    for ivar = 1:length(vars)
        cdata{ivar} = cat(2,cdata{ivar}, mydata{ivar});
    end
    cdata{ivar+1} = cat(2,cdata{ivar+1},files(ilog));
end

for ivar = 1:length(vars)
    vardata.(vars{ivar}) = cdata{ivar};
end
outdata.data = vardata;
outdata.logname = cdata{ivar+1};
outdata.vartype = vartype;
outdata.cfg = cfg;
