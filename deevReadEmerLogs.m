function outdata = deevReadEmerLogs(cfg)
% DeEv analysis script that reads in deev_*.txt files in 'logs/' from base
% emergent proj dir
%
%   input: 
%       cfg
%        filefiltstr: regexp filter on file listing, default is 'deev_36events_Sub*.txt'
%        badsubs: subject numbers to exclude from analysis
%        dir: struct containing the logs (def: logfiles)
%
%   output:
%       outdata: struc with headernames as matricies with Subs on the
%       columns



%defaults
if ~exist('cfg','var')          cfg = [];                end
if ~isfield(cfg,'badsubs')      cfg.badsubs = [];        end %which subjects to remove
if ~isfield(cfg,'filefiltstr')  cfg.filefiltstr = 'deev_36events_Sub[0-9]{1,2}\.txt';     end
if ~isfield(cfg,'dir')          cfg.dir = 'logs';        end %where are the log files?


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

%make sure logs are sorted by subject number
[goodsubs,ind] = sort(goodsubs);
files = files(ind);



%get var names and types from header
temp = textread(files{1},'%s');
vars = strread(temp{1},'%s','delimiter',',');
fstring = [];
vartype = {};
for ivar = 1:length(vars)
    switch vars{ivar}(1)
        case '_'
            %header
            fstring = [fstring '%s '];
            vartype{ivar} = '_';
            vars{ivar} = vars{ivar}(1:end-1);%remove colon
        case '$'
            fstring = [fstring '%s '];
            %outdata.(vars{ivar}(1:end-1)) = {};
            vartype{ivar} = '$';
        case '|'
            fstring = [fstring '%d '];
            %outdata.(vars{ivar}(1:end-1)) = [];
            vartype{ivar} = '|';
        case '%'
            fstring = [fstring '%f '];
            %outdata.(vars{ivar}(1:end-1)) = [];
            vartype{ivar} = '%';
        otherwise 
            fstring = [fstring '%s '];
            %outdata.(vars{ivar}(1:end-1)) = {};
            vartype{ivar} = '$';
    end
    vars{ivar} = vars{ivar}(2:end);    
end

%find c[0-9]t[0-9] vars
matind = regexp(vars,'^c[0-9]t[0-9]$');
matind = ~cellfun(@isempty,matind);

cdata = cell(1,length(vars)+1);
res = [];
for ilog = 1:length(files)    
    %read in log file
    logfilename = files{ilog};
    fid = fopen(logfilename,'r');
    mydata = textscan(fid, fstring, 'Headerlines', 1, 'Delimiter', ',', 'TreatAsEmpty', {'na'});        
    fclose(fid);
    
    catmat = [];
    for ivar = 1:length(vars)
        cdata{ivar} = cat(2,cdata{ivar}, mydata{ivar});
        if matind(ivar)             
            tmp = single(mydata{ivar});
            tmp(tmp==-1) = nan;
            catmat = cat(2,catmat,tmp);
        end
    end
    res = cat(3,res,catmat);
    cdata{ivar+1} = cat(2,cdata{ivar+1},files(ilog));
end

for ivar = 1:length(vars)
    vardata.(vars{ivar}) = cdata{ivar};
end
nEvts = size(res,1);
evtstr = strcat(repmat({'evt'},1,nEvts),strtrim(cellstr(num2str([1:nEvts]'))'));
substr = strcat(repmat({'sub'},1,length(goodsubs)),strtrim(cellstr(num2str(goodsubs'))'));

resDO{1} = evtstr;
resDO{2} = vars(matind');
resDO{3} = substr;

outdata.res = res;
outdata.resDO = resDO;
outdata.data = vardata;
outdata.lognames = cdata{ivar+1};
outdata.vartype = vartype;
outdata.cfg = cfg;
