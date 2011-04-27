function out = get(t,varargin)
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

if isa(t,'chanhandle')
    fields = fieldnames(t.handle(1));       % get list of fields
    switch nargin
    case 1        % display of all properties of object 
        for i = 1:length(t.handle)
            disp([ '%%% Properties of chanhandle #' num2str(i)]);
            for j = 1:length(fields)
                dispfield(t.handle(i),fields{j});
            end;
        end;
    case 2      % Return the value of only the specified property
        
        if isnumeric(varargin{1})   % Ex. get(t(channel_nums))
                                    % argin = channel number
            if min(varargin{1})>0 & ... % check if channel number is valid
                    max(varargin{1})<=length(t.handle) &...
                    (round(varargin{1})==varargin{1})
                chlist = varargin{1};
                out = chanhandle(t.handle(chlist(1)));
                for chi = 2:length(chlist)
                    out = [out;chanhandle(t.handle(chlist(chi)))];
                end
            else
                error('Channel number out of bounds');
            end;
        elseif ischar(varargin{1})  % Ex. get(t,'field')
                                    % argin = 'field' 
            prop = varargin{1};           
            if strcmpi(prop,'NumChannels') % If Prop = NumChannels
                out = length(t.handle);
                return;
            end;
            if ischar(prop),    % Check to see if the Property is valid
                testval = getsfield(t.handle(1),prop);
                if isnumeric(testval)   % If O/p value is numeric then simple array
                    for i = 1:length(t.handle)
                        out(i,:) = getsfield(t.handle(i),prop);
                    end;
                elseif ischar(testval)  % Elseif O/p value is string than cell array
                    for i = 1:length(t.handle)
                        out{i} = getsfield(t.handle(i),prop);
                    end;
                end;
            end;
        end;
    case 3  %Ex. get(t,channel_nums,field)
        if min(varargin{1})>0 & ...
                max(varargin{1})<=length(t.handle) &...
                (round(varargin{1})==varargin{1})       
            chnum = varargin{1};        
        else
            error('Channel number out of bounds');
        end;
        prop = varargin{2};
        if ischar(prop),
            testval = getsfield(t.handle(chnum(1)),prop);
            if isnumeric(testval)
                for i = 1:length(chnum)
                    out(i,:) = getsfield(t.handle(chnum(i)),prop);
                end;
            elseif ischar(testval)
                for i = 1:length(chnum)
                    out{i} = getsfield(t.handle(chnum(i)),prop);
                end;
            end;
        end;
    end;
end;
    
function stdout = dispfield(t,prop)
% DISPFIELD gets field value depending on field. Takes care of special
% situtations like for private fields, sensitivity and amplifier settings
% and acquisition parameters
% val = dispfield(t,prop)

if strcmpi(prop,'Sensitivity')
    % Sensitivity field is not transparent ot user. The properties
    % of Sensitivity are displayed as Senstivity.field
    stdout = sprintf('%s : %3.3f\n%s : %3.3f',...
        'Sensitivity.Gains',t.Sensitivity.Gains,...
        'Sensitivity.Offsets',t.Sensitivity.Offsets);    
elseif strcmpi(prop,'SensorInfo')
    % SensorInfo field is not transparent ot user. The properties
    % of SensorInfo are displayed as SensorInfo.field
    if max(strcmpi({'trigger','EEG','ECG','etc','NULL'},t.SensorInfo.type))
        stdout = sprintf('%s : %s\n%s : %s\n%s : %d',...
            'SensorInfo.Type',t.SensorInfo.type,...
            'SensorInfo.ChanName',t.SensorInfo.channame,...
            'SensorInfo.IdentNum',t.SensorInfo.identnum);
    else
        taglist = fieldnames(t.SensorInfo);
        stdout = sprintf('%s : %s\n','SensorInfo.Type',t.SensorInfo.type);
        taglist(find(strcmpi(taglist(:),'type'))) = [];
        for i = 1:length(taglist)
            stdout = [stdout,sprintf('%s : %3.3f\n',['SensorInfo.' taglist{i}],...
                    eval(['t.SensorInfo.' taglist{i}]))];
        end;
    end;
elseif strcmpi(prop,'Private')  % If private, display nothing
    return;
else                            % For all other properties
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
% situtations like for Sensitivity and SensorInfo
% val = getfield(t.handle,prop)

out = [];                   % init out
fields = fieldnames(t);     % Fields of substructure "handle"

subsindx = findstr(prop,'.');
if ~isempty(subsindx)               % if subfield of Sensitivity/SensorInfo
    subtag = prop(1:subsindx-1);    % Get main field (Sensitivity/SensorInfo)
    subsubs = prop(subsindx+1:end); % Get subfield
    ntag = find(strcmpi(fields(:),subtag));
    if ~isempty(ntag),              % If valid mainfield
        subtagstruct = getfield(t,fields{ntag}); % get structure of mainfield
        substags = fieldnames(subtagstruct);            
        nsubs = find(strcmpi(substags(:),subsubs));
        if ~isempty(nsubs),         % If valid subfield
            out = getfield(subtagstruct,substags{nsubs});
        end;
    end;
else                                % If main field
    n = find(strcmpi(fields(:),prop));
    if ~isempty(n), 
        out = getfield(t,fields{n});
    end;
end;  
if isempty(out)
    error('sorry - couldnt obtain channel info - please check input parameters');
end;



