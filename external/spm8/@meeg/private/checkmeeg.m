function [result meegstruct]=checkmeeg(meegstruct, option)
% Function for checking the internal struct of meeg objects
% FORMAT [result meegstruct]=checkmeeg(meegstruct, option)
% result - 1 - OK, 0- failed
% meegstruct - the struct to check (is returned modified if necessary)
% option - 'basic' (default) - just check the essential fields
%          'sensfid' - also checks sensor and fiducial definitions
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: checkmeeg.m 5026 2012-10-31 14:55:47Z vladimir $

if nargin==1
    option = 'basic';
end

result=0;

if ~isfield(meegstruct, 'Nsamples')
    disp('checkmeeg: number of samples per trial is missing');
    return;
else
    Nsamples = meegstruct.Nsamples;
end

if ~isfield(meegstruct, 'Fsample') && (Nsamples~=0)
    disp('checkmeeg: sampling rate is missing');
    return;
elseif (Nsamples==0)  % This is to enable creation of empty meeg objects
    meegstruct.Fsample = 0;
end

if ~isfield(meegstruct, 'timeOnset')
    meegstruct.timeOnset = 0;
elseif isempty(meegstruct.timeOnset)
    if meegstruct.Fsample ~= 0
        meegstruct.timeOnset = 1/meegstruct.Fsample;
    else
        meegstruct.timeOnset = 0; % This is to enable creation of empty meeg objects
    end
end


if ~isfield(meegstruct, 'trials') && (Nsamples~=0)
    disp('checkmeeg: no trials description');
    return;
elseif (Nsamples==0)  % This is to enable creation of empty meeg objects
    meegstruct.trials = struct([]);
else
    Ntrials = length(meegstruct.trials);
    if ~isfield(meegstruct.trials, 'label')
        disp('checkmeeg: no trial label, assigning default');
        [meegstruct.trials.label] = deal('Undefined');
    elseif any(cellfun('isempty', {meegstruct.trials(:).label}))
        disp('checkmeeg: some trial labels empty, assigning default');
        [meegstruct.trials(cellfun('isempty', {meegstruct.trials(:).label})).label] = deal('Undefined');
    end
    if ~isfield(meegstruct.trials, 'bad')
        [meegstruct.trials.bad] = deal(0);
    end
    if ~isfield(meegstruct.trials, 'events')
        [meegstruct.trials.events] = deal([]);
    end
    
    for k = 1:Ntrials
        
        label = meegstruct.trials(k).label;
        
        if iscell(label) && numel(label) == 1
            label = label{1};
        end
        
        if isnumeric(label)
            label = num2str(label);
        end
        
        if isa(label, 'char')
            meegstruct.trials(k).label = deblank(label);
        else
            meegstruct.trials(k).label = 'Unknown';
            disp('checkmeeg: some trial labels were not strings, changing back to ''Unknown''');
        end
        
        if  length(meegstruct.trials(k).bad)>1 || ~ismember(meegstruct.trials(k).bad, [0, 1])
            disp(['checkmeeg: illegal value for bad flag in trial ' num2str(k) ', setting to zero.']);
            meegstruct.trials(k).bad = 0;
        end
        
        event = meegstruct.trials(k).events;
        
        if ~isempty(event) && ~(numel(event) == 1 && isequal(event.type, 'no events'))
            % make sure that all required elements are present
            if ~isfield(event, 'type'),     error('type field not defined for each event');  end
            if ~isfield(event, 'time'),     error('time field not defined for each event');  end
            if ~isfield(event, 'value'),    [event.value]    = deal([]);                     end
            if ~isfield(event, 'offset'),   [event.offset]   = deal(0);                      end
            if ~isfield(event, 'duration'), [event.duration] = deal([]);                     end
            
            
            % make sure that all numeric values are double
            for i=1:length(event)
                if isnumeric(event(i).value)
                    event(i).value = double(event(i).value);
                end
                event(i).time      = double(event(i).time);
                event(i).offset    = double(event(i).offset);
                event(i).duration  = double(event(i).duration);
            end
            
            if ~isempty(event)
                % sort the events on the sample on which they occur
                % this has the side effect that events without time are discarded
                [dum, indx] = sort([event.time]);
                event = event(indx);
            end
        else
            event = [];
        end
        
        meegstruct.trials(k).events = event;
    end
    
    if ~isfield(meegstruct.trials, 'onset')
        [meegstruct.trials.onset] = deal(0);
    else
        [meegstruct.trials(find(cellfun('isempty', {meegstruct.trials.onset}))).onset] = deal(0);
    end
    if ~isfield(meegstruct.trials, 'repl')
        [meegstruct.trials.repl] = deal(1);
    end
end

if ~isfield(meegstruct, 'channels') && (Nsamples~=0)
    disp('checkmeeg: no channels description');
    return;
elseif (Nsamples==0)  % This is to enable creation of empty meeg objects
    meegstruct.channels = struct([]);
else
    Nchannels = length(meegstruct.channels);
    if ~isfield(meegstruct.channels, 'label')
        disp('checkmeeg: no channel label, assigning default');
        for i = 1:Nchannels
            meegstruct.channels(i).label = ['Ch' num2str(i)];
        end
    end
    if ~isfield(meegstruct.channels, 'bad')
        [meegstruct.channels.bad] = deal(0);
    else
        [meegstruct.channels(find(cellfun('isempty', {meegstruct.channels.bad}))).bad] = deal(0);
    end
    
    for i = 1:Nchannels
        meegstruct.channels(i).bad = double(~~meegstruct.channels(i).bad);
    end
    
    if ~isfield(meegstruct.channels, 'type')
        disp('checkmeeg: no channel type, assigning default');
        [meegstruct.channels.type] = deal('Other');
    end
    
    % This is for backward compatibility with early SPM8b. Can be removed
    % after a while.
    ind = strmatch('MEGREF', {meegstruct.channels.type}, 'exact');
    [meegstruct.channels(ind).type] = deal('REF');
    
    if ~isfield(meegstruct.channels, 'X_plot2D')
        [meegstruct.channels.X_plot2D] = deal([]);
        [meegstruct.channels.Y_plot2D] = deal([]);
    end
    if ~isfield(meegstruct.channels, 'units')
        disp('checkmeeg: no units, assigning default');
        [meegstruct.channels.units] = deal('unknown');
    else
        [meegstruct.channels(find(cellfun('isempty', {meegstruct.channels.units}))).units] = deal('unknown');
    end
end

try
    meegstruct.transform.ID;
catch
    meegstruct.transform.ID = 'time';
    disp('checkmeeg: transform type missing, assigning default');
end

if strncmp(meegstruct.transform.ID, 'TF', 2) % TF and TFphase
    try
        Nfrequencies = length(meegstruct.transform.frequencies);
    catch
        error('Information about frequencies missing');
    end
end

if isfield(meegstruct, 'path')
    filepath = meegstruct.path;
else
    filepath = pwd;
end

if ~isfield(meegstruct, 'data') && (Nsamples~=0)
    disp('checkmeeg: no data field');
    return;
elseif (Nsamples==0)  % This is to enable creation of empty meeg objects
    meegstruct.data = struct([]);
else
    if ~isfield(meegstruct.data, 'fnamedat')
        if ~isa(meegstruct.data, 'file_array')
            disp('checkmeeg: data file name missing');
            return;
        else
            disp('checkmeeg: SPM12 format - converting back');
            sav_sc = meegstruct.data.scl_slope;
            sav_os = meegstruct.data.offset;
            meegstruct.data = struct('fnamedat', meegstruct.data.fname, ...
                'datatype', meegstruct.data.dtype);
        end
    end
    
    fname = meegstruct.data.fnamedat;
    if ~ispc, fname = strrep(fname,'\',filesep); end
    [junk, fnamedat, ext] = fileparts(fname);
    if isempty(ext)
        meegstruct.data.fnamedat = [fnamedat '.dat'];
    else
        meegstruct.data.fnamedat = [fnamedat ext];
    end
    
    if ~isfield(meegstruct.data, 'datatype')
        disp('checkmeeg: data type missing, assigning default');
        meegstruct.data.datatype = 'float32';
    end
    
    if ~isfield(meegstruct.data, 'y')
        meegstruct.data.y=[];
    end
    
    if isa(meegstruct.data.y, 'file_array')
        % catching up (unlikely case) where filearray.fname is
        % different from data.fnamedat -> set data.fnamedat
        fname = meegstruct.data.y.fname;
        if ~ispc, fname = strrep(fname,'\',filesep); end
        [junk, yfname, yext] = fileparts(fname);
        
        fnamedat = meegstruct.data.fnamedat;
        if ~ispc, fnamedat = strrep(fnamedat,'\',filesep); end
        [junk, dfname, dext] = fileparts(fnamedat);
        if ~strcmp([yfname yext],[dfname dext])
            meegstruct.data.fnamedat = [yfname yext];
        end
        
        meegstruct.data.y.fname = fullfile(filepath, [yfname yext]);
        
        % save original file_array scale & offset, just in case
        sav_sc = meegstruct.data.y.scl_slope;
        sav_os = meegstruct.data.y.offset;
        try
            % Try reading data, i.e. check if it's a "good" filearray
            meegstruct.data.y(1, 1, 1);
        catch
            meegstruct.data.y = [];
        end
    end
    
    if isa(meegstruct.data.y, 'file_array')
        fn = meegstruct.data.y.fname;
        if ~ispc, fn = strrep(fn,'\',filesep); end
    else
        fn = '';
    end
    if ~isa(meegstruct.data.y, 'file_array') ...
            || isempty(fileparts(fn)) ...
            || ~exist(meegstruct.data.y.fname, 'file')
        switch(meegstruct.transform.ID)
            case 'time'
                meegstruct.data.y = file_array(fullfile(filepath, meegstruct.data.fnamedat), ...
                    [Nchannels Nsamples Ntrials], meegstruct.data.datatype);
                
            case {'TF', 'TFphase'}
                meegstruct.data.y = file_array(fullfile(filepath, meegstruct.data.fnamedat), ...
                    [Nchannels Nfrequencies Nsamples Ntrials], meegstruct.data.datatype);
                if Ntrials>1
                    expected_size = [Nchannels Nfrequencies Nsamples Ntrials];
                else
                    expected_size = [Nchannels Nfrequencies Nsamples];
                end
                
            otherwise
                error('Unknown transform type');
        end
        % and restore original file_array scale, if available (exist) & useful (~=[])
        if exist('sav_sc','var') && ~isempty(sav_sc) && ...
                size(meegstruct.data.y, 1) == length(sav_sc)
            meegstruct.data.y.scl_slope = sav_sc;
        end
        % and restore original file_array offset, if available (exist) & useful (~=0)
        if exist('sav_os','var') && sav_os
            meegstruct.data.y.offset = sav_os;
        end
        
    end
    
    switch(meegstruct.transform.ID)
        case 'time'
            if Ntrials>1
                expected_size = [Nchannels Nsamples Ntrials];
            else
                expected_size = [Nchannels Nsamples];
            end
        case {'TF', 'TFphase'}
            if Ntrials>1
                expected_size = [Nchannels Nfrequencies Nsamples Ntrials];
            else
                expected_size = [Nchannels Nfrequencies Nsamples];
            end
            
            
        otherwise
            error('Unknown transform type');
    end
    
    if any(size(meegstruct.data.y) ~= expected_size)
        disp('checkmeeg: data size does not match the header');
        return;
    end
end

if ~isfield(meegstruct, 'type') ||...
        (strcmp(meegstruct.type, 'continuous') && Ntrials>1) ||...
        (strcmp(meegstruct.type, 'evoked') && (numel(unique({meegstruct.trials.label})) ~= Ntrials)) ||...
        (strcmp(meegstruct.type, 'continuous') && strncmp(meegstruct.transform.ID, 'TF', 2)) ||...
        strcmp(meegstruct.type, 'grandmean')
    disp('checkmeeg: data type is missing or incorrect, assigning default');
    % rule of thumb - 10 sec
    if Nsamples == 0
        meegstruct.type = 'continuous';
    elseif Ntrials==1 && (Nsamples/meegstruct.Fsample) > 10 &&...
            ~strncmp(meegstruct.transform.ID, 'TF', 2)
        meegstruct.type = 'continuous';
    elseif numel(unique({meegstruct.trials.label})) == Ntrials
        meegstruct.type = 'evoked';
    else
        meegstruct.type = 'single';
    end
end

try
    fname = meegstruct.data.y.fname;
    if ~ispc, fname = strrep(fname,'\',filesep); end
    [pdat, fdat] = fileparts(fname);
catch
    fdat = 'spm8';
end

if ~isfield(meegstruct, 'fname')
    meegstruct.fname = [fdat '.mat'];
else
    fname = meegstruct.fname;
    if ~ispc, fname = strrep(fname,'\',filesep); end
    [p, f] = fileparts(fname);
    if isempty(f)
        f = fdat;
    end
    meegstruct.fname = [f '.mat'];
end

if ~isfield(meegstruct, 'path')
    try
        fname = meegstruct.data.y.fname;
        if ~ispc, fname = strrep(fname,'\',filesep); end
        meegstruct.path = fileparts(fname);
    catch
        meegstruct.path = pwd;
    end
end

if ~isfield(meegstruct, 'sensors')
    meegstruct.sensors = struct([]);
else
    if isfield(meegstruct.sensors, 'eeg')
        if isempty(meegstruct.sensors.eeg)
            meegstruct.sensors = rmfield(meegstruct.sensors, 'eeg');
        else
            meegstruct.sensors.eeg = ft_datatype_sens(meegstruct.sensors.eeg);
        end
    end
    if isfield(meegstruct.sensors, 'meg')
        if isempty(meegstruct.sensors.meg)
            meegstruct.sensors = rmfield(meegstruct.sensors, 'meg');
        else
            meegstruct.sensors.meg = ft_datatype_sens(meegstruct.sensors.meg);
        end
    end
end

if ~isfield(meegstruct, 'fiducials')
    meegstruct.fiducials = struct([]);
end

if ~isfield(meegstruct, 'artifacts')
    meegstruct.artifacts = struct([]);
end

if ~isfield(meegstruct, 'transform')
    meegstruct.transform = struct('ID', 'time');
end

if ~isfield(meegstruct, 'other')
    meegstruct.other = struct([]);
else
    if isfield(meegstruct.other, 'origchantypes')
        % This is for backward compatibility
        % Can be removed after a while
        if iscell(meegstruct.other.origchantypes)
            if numel(meegstruct.other.origchantypes) == Nchannels
                origchantypes = meegstruct.other.origchantypes;
                meegstruct.other.origchantypes = struct([]);
                meegstruct.other.origchantypes(1).type = origchantypes;
                meegstruct.other.origchantypes(1).label = {meegstruct.channels(:).label}';
            else
                meegstruct.other = rmfield(meegstruct.other, 'origchantypes');
            end
        end
    end
    if isfield(meegstruct.other, 'condlist') && ~isempty(meegstruct.other.condlist)
        % This a workaround for the absence of 'first' option in unique()
        % of Matlab 7.1 and last occurence being the default
        clist = meegstruct.other.condlist(end:-1:1);
        [junk, ind] = unique(clist);
        clist = clist(sort(ind));
        meegstruct.other.condlist = clist(end:-1:1);
    end
end

if ~isfield(meegstruct, 'history')
    meegstruct.history = struct([]);
end

if ~isfield(meegstruct, 'cache')
    meegstruct.cache = struct([]);
end

% This makes sure the order of the fields in the struct is right. This
% seems to matter to the class function.

fieldnames_order = {
    'type'
    'Nsamples'
    'Fsample'
    'timeOnset'
    'trials'
    'channels'
    'data'
    'fname'
    'path'
    'sensors'
    'fiducials'
    'artifacts'
    'transform'
    'other'
    'history'
    'cache'};

[sel1, sel2] = match_str(fieldnames_order, fieldnames(meegstruct));
tempcell = struct2cell(meegstruct);
meegstruct = cell2struct(tempcell(sel2), fieldnames_order, 1);

if strcmp(option, 'basic')
    result = 1;
    return;
end

chantypes = getset(meegstruct, 'channels', 'type');
eegind = strmatch('EEG', chantypes, 'exact');
megind = strmatch('MEG', chantypes);
lfpind = strmatch('LFP', chantypes, 'exact');

% Allow DCM on a pure LFP dataset
if strcmp(option, 'dcm') && isempty([eegind(:); megind(:)])...
        && ~isempty(lfpind) && ismember(meegstruct.transform.ID, {'time', 'TF'})
    result = 1;
    return;
end

if strcmp(option, 'sensfid') || strcmp(option, '3d') ||...
        (strcmp(option, 'dcm') && ~isempty([eegind(:); megind(:)]))
    
    if ~strcmp(meegstruct.transform.ID, 'time')
        disp('checkmeeg: incorrect data type for source reconstruction');
        return;
    end
    
    if isempty(meegstruct.sensors)
        disp('checkmeeg: no sensor positions are defined');
        return;
    end
    
    if ~isempty(eegind)
        if ~isfield(meegstruct.sensors, 'eeg') || isempty(meegstruct.sensors.eeg)
            disp('checkmeeg: EEG channel locations are not specified');
            return;
        else
            if ~isempty(setdiff({meegstruct.channels(eegind).label}, meegstruct.sensors.eeg.label))
                disp('checkmeeg: not all EEG channel locations are specified');
                return;
            end
        end
    end
    
    if ~isempty(megind)
        if ~isfield(meegstruct.sensors, 'meg') || isempty(meegstruct.sensors.meg)
            disp('checkmeeg: MEG channel locations are not specified');
            return;
        else
            if ~isempty(setdiff({meegstruct.channels(megind).label}, meegstruct.sensors.meg.label))
                disp('checkmeeg: not all MEG channel locations are specified');
                return;
            end
        end
    end
    
    if isempty(meegstruct.fiducials)
        disp('checkmeeg: no fiducials are defined');
        return;
    end
    
    if ~isfield(meegstruct.fiducials, 'pnt') || isempty(meegstruct.fiducials.pnt)
        if ~isempty(eegind)
            % Copy EEG sensors to fiducials.
            meegstruct.fiducials.pnt = meegstruct.sensors.eeg.elecpos;
        else
            meegstruct.fiducials.pnt = sparse(0, 3);
        end
    end
    
    if ~isfield(meegstruct.fiducials, 'fid') || ...
            ~all(isfield(meegstruct.fiducials.fid, {'pnt', 'label'})) ||...
            (length(meegstruct.fiducials.fid.label) ~= size(meegstruct.fiducials.fid.pnt, 1)) || ...
            length(meegstruct.fiducials.fid.label) < 3
        disp('checkmeeg: at least 3 fiducials with labels are required');
        return
    end
    
    nzlbl = {'fidnz', 'nz', 'nas', 'nasion', 'spmnas'};
    lelbl = {'fidle', 'fidt9', 'lpa', 'lear', 'earl', 'le', 'l', 't9', 'spmlpa'};
    relbl = {'fidre', 'fidt10', 'rpa', 'rear', 'earr', 're', 'r', 't10', 'spmrpa'};
    
    [sel1, nzind] = match_str(nzlbl, lower(meegstruct.fiducials.fid.label));
    if isempty(nzind)
        disp('checkmeeg: could not find the nasion fiducial');
    else
        nzind = nzind(1);
    end
    
    [sel1, leind] = match_str(lelbl, lower(meegstruct.fiducials.fid.label));
    if isempty(leind)
        disp('checkmeeg: could not find the left fiducial');
    else
        leind = leind(1);
    end
    
    [sel1, reind] = match_str(relbl, lower(meegstruct.fiducials.fid.label));
    if isempty(reind)
        disp('checkmeeg: could not find the right fiducial');
    else
        reind = reind(1);
    end
    
    restind = setdiff(1:length(meegstruct.fiducials.fid.label), [nzind(:)', leind(:)', reind(:)']);
    
    meegstruct.fiducials.fid.label = meegstruct.fiducials.fid.label([nzind(:)', leind(:)', reind(:)', restind(:)']);
    meegstruct.fiducials.fid.pnt = meegstruct.fiducials.fid.pnt([nzind(:)', leind(:)', reind(:)', restind(:)'], :);
    
    result = 1;
end



