function out = get(t,prop)
% GET    Get object properties.
% out = GET(t,'PropertyName') returns the value of the specified
% property for the object t.  If t is a vector of handles, then 
% get will return an M-by-1 cell array of values where M is equal
% to length(t).  If 'PropertyName' is replaced by a 1-by-N or N-by-1
% cell array of strings containing property names, then GET will 
% return an M-by-N cell array of values.
%  
% GET(t) displays all property names and values

if nargin<1,
    error('first argument must be acqparam object')
end;

if isa(t,'sqdhandle')
    fields = fieldnames(t); % get fields of object
    switch nargin
    case 1  % Display all properties
        for j = 1:length(fields)
            dispfield(t,fields{j});
        end;
    case 2 % return specified property
        out = [];
        if ischar(prop),
            out = getsfield(t,prop);
        elseif iscell(prop)
            for j = 1:length(prop)
                out{j} = getsfield(t,prop{j});
            end;
        end;
        
    end;
end;

function stdout = dispfield(t,prop)
% DISPFIELD gets field value depending on field. Takes care of special
% situtations like for private fields, sensitivity and amplifier settings
% and acquisition parameters
% val = dispfield(t,prop)

if strcmpi(prop,'acqparam')
    get(t.acqparam);
    stdout = '';
elseif strcmpi(prop,'Channels')
    % For Channels, just display [1xChannelNumer] array of chanhandle objects
    stdout = sprintf('%s : %s','Channels',...
        ['[1x' num2str(get(t,'ChannelCount')) '] Channel array']);
elseif strcmpi(prop,'PatientInfo')
    % For patient info, just display Name of the patient person
    stdout     = sprintf('%s : %s','PatientInfo',...
        t.PatientInfo.name);   
elseif strcmpi(prop,'Amplifier')
    % The field Amplifier is transparent ot the user. Directly display
    % the inputgain and output gain
    stdout     = sprintf('%s : %d\n%s : %d',...
        'InputGain',t.Amplifier.InputGain,...
        'OutputGain',t.Amplifier.OutputGain);
elseif strcmpi(prop,'FileName')
    % Display filename without path
    [pathstr,fname] = fileparts(t.FileName);
    stdout = sprintf('%s : %s','FileName',fname);
elseif strcmpi(prop,'Private')  % if private
    return;                     % show nothing
else
    % For all otherproperties
    val = getsfield(t,prop);
    val = val';
    if ischar(val),
        stdout = sprintf('%s : %s',prop,val);
    elseif max(val-floor(val))>0
        stdout = sprintf('%s : %s',prop,sprintf('%10.3f',val(:)));
    else
        stdout = sprintf('%s : %s',prop,sprintf('%d',val(:)));
    end;
end;
disp(stdout);

function out = getsfield(t,prop)
% GETSFIELD gets field value depending on field. Takes care of special
% situtations like for patientinfo,amplifier settings, 
% and acquisition parameters
% val = getfield(t,prop)
out = [];
taglist = fieldnames(struct(t));  % List of fields  
if strcmpi(prop,'InputGain')        % if Amplifier.InputGain
    out = t.Amplifier.InputGain;
elseif strcmpi(prop,'OutputGain')   % if Amplifier.OutputGain
    out = t.Amplifier.OutputGain;
elseif strcmpi(prop,'SensitivityGain')
    sensinfo = sensitivityinfo(t.FileName);
    out = sensinfo.Gains;
else                                % for others
    n = find(strcmpi(taglist(:),prop));
    if ~isempty(n),                 % If main field
        out = getfield(struct(t),taglist{n});
    else                            % if patient info
        patinfolist = fieldnames(t.PatientInfo);
        subindx = find(strcmpi(patinfolist(:),prop));
        if ~subindx
            out = getfield(t.PatientInfo,patinfolist{subindx});
        else                        % if property of acqparam
            out = get(t.acqparam,prop);
        end;
    end;
end;




