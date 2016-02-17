function trl = deev_trialfun_mff(cfg)

% operates using Net Station evt files and event structs

% convert single string into cell-array, otherwise intersection does not
% work as intended
if ischar(cfg.trialdef.eventvalue)
    cfg.trialdef.eventvalue = {cfg.trialdef.eventvalue};
end


% get the header and event information
fprintf('Reading flags from EEG file using FieldTrip...');
ft_hdr = ft_read_header(cfg.dataset);
[pathstr,name] = fileparts(cfg.dataset);
ftEventsFile = fullfile(pathstr,sprintf('%s_ftEvents.mat',name));
if exist(ftEventsFile,'file')
    ft_event = load(ftEventsFile);
    if isfield(ft_event,'date_string')
        warning('Using pre-saved FT events from this date: %s!',ft_event.date_string);
    else
        warning('Using pre-saved FT events from an unknown date!');
    end
    ft_event = ft_event.ft_event;
else
    tic
    ft_event = ft_read_event(cfg.dataset);
    toc
    date_string = datestr(now);
    fprintf('Saving FT events from MFF (current time: %s): %s...',date_string,ftEventsFile);
    save(ftEventsFile,'ft_event','date_string');
    fprintf('Done.\n');
end
fprintf('Done.\n');

%convert event types that start with ECI to be name 'TCPIP'
%hack to fix bad characters in netstation 5.2 ECI track name
types = {ft_event.type};
ind = find(~cellfun(@isempty,regexp(types,'ECI.*')));
for i = ind
    ft_event(i).type = 'TCPIP';
end
cfg.trialdef.eventtype = {'TCPIP'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in external data, if wanted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg.eventinfo.useMetadata
    md = cfg.eventinfo.metadata;
    
    if ismember('eventStruct', md.types)
        % read the events file
        eventsFile = fullfile(md.dirs.dataroot,md.dirs.behDir,md.subject,'events','events.mat');
        if exist(eventsFile,'file')
            fprintf('Loading events file: %s...',eventsFile);
            events_all = load(eventsFile,'events');
            events_all = events_all.events;
            fprintf('Done.\n');
        else
            error('Cannot find events file: %s\n',eventsFile)
        end
    end
    
    if ismember('behData', md.types)
        %read behavioral data for this subject
        cfg_beh = [];
        cfg_beh.dir = fullfile(md.dirs.dataroot, md.dirs.behDir);
        subnum = regexp(md.subject,'_','split'); subnum = subnum{end};
        cfg_beh.filefiltstr = ['deevEEGExp_1EEG_4blocks_36events_Sub' subnum '.txt'];
        behdata = deevBlkGetLogDep(cfg_beh);
        %get average accuracy across retrievals
        avgacc = squeeze(nanmean(behdata.res,2));
    end
    
    if ismember('expParam', md.types)
        % read the experiment parameters file
        expParamFile = fullfile(md.dirs.dataroot,md.dirs.behDir,md.subject,'experimentParams.mat');
        if exist(expParamFile,'file')
            fprintf('Loading experiment parameters file: %s...',expParamFile);
            load(expParamFile,'expParam');
            fprintf('Done.\n');
        else
            error('Cannot find experiment parameters file: %s\n',expParamFile)
        end
    end
end

%define triggers to look into
triggers = {'WRD1','TRG2'};
if cfg.eventinfo.usePhotodiodeDIN
    triggers = cat(2,triggers,cfg.eventinfo.photodiodeDIN_str);
    %cfg.trialdef.eventtype = cat(2, cfg.trialdef.eventtype, cfg.eventinfo.photodiodeDIN_str);
end

% initialize the trl matrix
trl = [];

% all trls need to have the same length
maxTrlCols = -Inf;
fn_trl_ord = fieldnames(cfg.eventinfo.trl_order);
for fn = 1:length(fn_trl_ord)
    if ismember(fn_trl_ord{fn},cfg.eventinfo.eventValues)
        if length(cfg.eventinfo.trl_order.(fn_trl_ord{fn})) > maxTrlCols
            maxTrlCols = length(cfg.eventinfo.trl_order.(fn_trl_ord{fn}));
        end
    end
end
if maxTrlCols == -Inf
    error('Did not set maximum number of trialinfo columns!\n');
end
timeCols = 3;
eventNumCols = 1;
trl_ini = -1 * ones(1, timeCols + eventNumCols + maxTrlCols);

if cfg.eventinfo.usePhotodiodeDIN
    photodiodeDIN_toleranceMS = cfg.eventinfo.photodiodeDIN_toleranceMS;
    photodiodeDIN_toleranceSamp = ceil((photodiodeDIN_toleranceMS / 1000) * ft_hdr.Fs);        
end

if ~isfield(cfg.eventinfo,'offsetMSTCP')
    error('need to set default TCP offset!');
else
    offsetTCPSamp = ceil((cfg.eventinfo.offsetMSTCP/1000)*ft_hdr.Fs);
end

if cfg.eventinfo.offsetMS == -1
    %auto set anti-aliasing offset based on (note: assuming NA400 amp)
    % http://hosted.verticalresponse.com/662876/25e92832cf/1469694273/6774d3b846/
    cfg.eventinfo.offsetMS = 0;
    if str2num(ft_hdr.orig.xml.info.acquisitionVersion(1:3))<5.1 %5.2 corrects for this offset automatically
        switch ft_hdr.Fs
            case 1000
                cfg.eventinfo.offsetMS =  36;
            case 500
                cfg.eventinfo.offsetMS = 66;
            case 250
                cfg.eventinfo.offsetMS = 112;
        end
    end
end
offsetSamp = ceil((cfg.eventinfo.offsetMS / 1000) * ft_hdr.Fs);

%remove types with missing values
ft_event = ft_event(~cellfun(@isempty,{ft_event.value}));

% only keep the ft events with triggers
ft_event = ft_event(ismember({ft_event.value},triggers));

%any DIN triggers left?
if cfg.eventinfo.usePhotodiodeDIN
    isDIN = sum(strcmp(cfg.eventinfo.photodiodeDIN_str,{ft_event.value}))>0;
    if ~isDIN
        warning('No photodiode events found with value %s\ndefault TCP offset of %.2fms will be used',...
            cfg.eventinfo.photodiodeDIN_str, cfg.eventinfo.offsetMSTCP);
    end
end

%% go through events and add metadata to trl matrix

ses = cfg.eventinfo.sessionNum;
sesName = cfg.eventinfo.sessionNames{ses};
sesType = find(ismember(cfg.eventinfo.sessionNames,cfg.eventinfo.sessionNames{ses}));
% sesType = ismember(cfg.eventinfo.sessionNames,cfg.eventinfo.sessionNames{ses});

fprintf('FT event count of NS flags (out of %d): %s',length(ft_event),repmat(' ',1,length(num2str(length(ft_event)))));

for i = 1:length(ft_event)
    fprintf(1,[repmat('\b',1,length(num2str(i))),'%d'],i);
    
    if strcmp(ft_event(i).type,cfg.trialdef.eventtype)
        % found an EEG event that we might want to process
        nKeys = length(ft_event(i).orig.keys);
        
        % set column types because Net Station evt files can vary
        ns_evt_cols = {};
        for ns = 1:nKeys
            ns_evt_cols = cat(1,ns_evt_cols,ft_event(i).orig.keys(ns).key.keyCode);
        end
        
        %find which stim categories are presented
        cols.ordn = find(strcmp(ns_evt_cols,'ordn'));
        ordn = regexp(ft_event(i).orig.keys(cols.ordn).key.data.data,'-','split');
        
        keeptrial = 0;
        switch ft_event(i).value
            case 'WRD1'
                keeptrial = 1;
                %encoding trials
                
                %coded stim types (i.e. 1 through 4) to match ordn coding
                stimtypes = {'encLoc','encPer','encObj','encAni'}; 
                notstimtypes = {'encNotLoc','encNotPer','encNotObj','encNotAni'};
                evVal = stimtypes(cellfun(@str2num,ordn));
                evVal = cat(2,evVal,notstimtypes(~ismember(stimtypes,evVal)));

            case 'TRG2'
                keeptrial = 1;
                %retrieval trials
                
                %is this OL or CL?
                cols.cond = find(strcmp(ns_evt_cols,'sttp'));
                if isempty(cols.cond)
                    error('could not find stim type column to determine condition')
                end
                
                sttp = ft_event(i).orig.keys(cols.cond).key.data.data;
                evVal = {sprintf(sttp)};
                                
        end %switch
        
        if keeptrial
            
            % find where this event type occurs in the list
            eventNumber = find(ismember(cfg.trialdef.eventvalue,evVal));
            if isempty(eventNumber)
                error('event number not found for %s!\n',evVal{1});
            end
            
            fulleventNumber = eventNumber;
            fullevVal = evVal;
            evtcnt = 0;
            for eventNumber = fulleventNumber
                evtcnt = evtcnt + 1;
                evVal = fullevVal{evtcnt}; 
                                                
                % set the times we need to segment before and after the
                % trigger
                prestimSec = abs(cfg.eventinfo.prepost{1}(eventNumber,1));
                poststimSec = cfg.eventinfo.prepost{1}(eventNumber,2);
                
                % prestimulus period should be negative
                prestimSamp = -round(prestimSec * ft_hdr.Fs);
                poststimSamp = round(poststimSec * ft_hdr.Fs);
                
                
                % add it to the trial definition
                this_trl = trl_ini;
                
                % get the time of this event
                this_sample = ft_event(i).sample;
                
                % if we're using the photodiode DIN and we find one
                % within the threshold, replace the current sample time
                % with that of the DIN
                % otherwise use default offset calculated from timing tests
                if cfg.eventinfo.usePhotodiodeDIN && isDIN && (i~=1 && i<length(ft_event))
                    
                    %find time between next and previous event
                    diffs = [ft_event(i+1).sample ft_event(i-1).sample] - this_sample;
                    %are they DIN events?
                    DINstr = ismember({ft_event(i+1).value ft_event(i-1).value},cfg.eventinfo.photodiodeDIN_str);
                    %set non-DIN events to infinity
                    diffs(~DINstr) = Inf;
                    %take DIN event with minimum time from TCP event
                    [~,imin] = min(abs(diffs));
                    
                    mydiff = diffs(imin);
                    myind = [1 -1]; myind = myind(imin);
                    
                    %is it in our acceptable range
                    if abs(mydiff) <= photodiodeDIN_toleranceSamp
                        this_sample = ft_event(i+myind).sample;
                    else
                        %min DIN not in tolerance range
                        %use default offset here
                        warning('min DIN not in tolerance range\nusing default TCP offset of %.2fms',cfg.eventinfo.offsetMSTCP);
                        this_sample = ft_event(i).sample + offsetTCPSamp;
                    end
                else
                    %use default tcp offset
                    this_sample = ft_event(i).sample + offsetTCPSamp;                    
                end
                
                % prestimulus sample
                this_trl(1) = this_sample + prestimSamp + offsetSamp;
                % poststimulus sample
                this_trl(2) = this_sample + poststimSamp + offsetSamp;
                % offset in samples
                this_trl(3) = prestimSamp;
                
                %event number
                this_trl(4) = eventNumber;
                
                %add trial info
                % 'subn', 'exbk', 'phsn', 'trln', 'evtn', 'blkn', 'ordn',...
                %   'sttp', 'pstn', 'ltwn', 'ltwn', 'rtmm', 'rtcf', 'conf', 'crct';
                trl_order = cfg.eventinfo.trl_order.(evVal);
                for to = 1:length(trl_order)
                    thisInd = find(ismember(trl_order,trl_order{to}));
                    keyInd = find(ismember(ns_evt_cols, trl_order{to}));
                    
                    if ~isempty(thisInd) & ~isempty(keyInd)
                        evtvar = [];
                        switch trl_order{to}
                            case {'subn', 'exbk', 'phsn', 'trln', 'evtn', 'blkn',...
                                    'pstn', 'ltwn', 'rtmm', 'rtcf', 'conf', 'crct'}
                                %numbers
                                evtvar = str2num(ft_event(i).orig.keys(keyInd).key.data.data);
                                
                            case {'ordn'} %numbers with hyphen, just remove hypen
                                evtvar = str2num(cell2mat(ordn));
                                
                            case {'sttp'}
                                %strings
                                stimtypes = {'OL','CL'};%1=openloop, 2=closedloop
                                evtvar = find(ismember(stimtypes,ft_event(i).orig.keys(keyInd).key.data.data));
                        end
                        
                        this_trl(timeCols + eventNumCols + thisInd) = evtvar;
                    end
                    
                end %for each trial order
                
                % put all the trials together
                trl = cat(1,trl,double(this_trl));
            end %for each event number
        end % if keeptrial
    end
end

fprintf('\n');
