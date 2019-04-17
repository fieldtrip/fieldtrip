function out = get(t,prop)
% GET(T,PROP)    Get object properties.
% out = GET(t,'PropertyName') returns the value of the specified
% property for the object t.  If t is a vector of handles, then 
% get will return an M-by-1 cell array of values where M is equal
% to length(t).  If 'PropertyName' is replaced by a 1-by-N or N-by-1
% cell array of strings containing property names, then GET will 
% return an M-by-N cell array of values.
%  
% GET(t) displays all property names and values

if nargin<1,
    error('first argument must be acqparam object');
end;

if ~isa(t,'acqparam')
    error('first argument must be acqparam object');
else
    taglist = fieldnames(struct(t));    % Get list of fields of object
    switch nargin
    case 1      % If just asking for a regular display of all properties
        for j = 1:length(taglist)
            dispfield(t,taglist{j});
        end;
    case 2      % If asking to return value of specific properties
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

% To Display SampleInfo:
% The field SampleInfo is transparent to the user
% The user is shown directly the contents of SampleInfo
sampinfotaglist = fieldnames(t.SampleInfo); % get list of fields of SampleINfo
if strcmpi(prop,'SampleInfo')               % To display whole of SampleInfo
    stdout = '';                            % Init output
    for i = 1:length(sampinfotaglist)
        stdout = [stdout,...                % show fields of SampleInfo
                sprintf('%s : %d\n',sampinfotaglist{i},...
                getfield(t.SampleInfo,sampinfotaglist{i}))];
    end;
    
elseif find(strcmpi(sampinfotaglist(:),prop))% To display specific property
                                             % of SampleInfo    
   indx = find(strcmpi(sampinfotaglist(:),prop)); 
   stdout = sprintf('%s : %d\n',sampinfotaglist{indx},...
                getfield(t.SampleInfo,sampinfotaglist{indx}));
            
elseif strcmpi(prop,'Private')              % If private property
    stdout = '';                            % show nothing
    
else                                        % For all other properties
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
% situtations like for SampleInfo
% val = getsfield(t,prop)
taglist = fieldnames(struct(t));            % Get regular fields of object
sampinfotaglist = fieldnames(t.SampleInfo); % Get fields of SampleInfo

if find(strcmpi(taglist(:),prop))           % If regular field
    indx = find(strcmpi(taglist(:),prop));
    out = getfield(struct(t),taglist{indx}); % see builtin getfield for structure
elseif find(strcmpi(sampinfotaglist(:),prop))% Is SampleInfo field
    indx = find(strcmpi(sampinfotaglist(:),prop));
    out = getfield(t.SampleInfo,sampinfotaglist{indx});
end;
