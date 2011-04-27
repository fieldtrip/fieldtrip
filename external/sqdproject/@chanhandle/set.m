function t = set(t,varargin)
% SET Set object properties
% Usage:
% t = set(t,property1,value1,property2,value2,...)
%	Sets values of properties 'property1', 'property2', etc. for chanhandle object 't'
%   to 'value1', value2', ... respectively 

if nargin<1,
    error('first argument must be chanhandle object')
end;

if ~isa(t,'chanhandle')
    error('First input must belong to class chanhandle');
end;
if ~rem(nargin-1,2)             % Ex. t = set(t,'prop',val);
    inargs = varargin;          % prop and val must come in pairs
    taglist = fieldnames(t.handle(1));
    chnum = 1:length(t.handle); 
else    
    error('Property and values must be given as pairs');
end;

if size(chnum,2)==1
    chnum = chnum'; % Need row vector
end;
for i = chnum
    while length(inargs)>=2
        prop = inargs{1};   % one property at a time
        val  = inargs{2};
        inargs = inargs(3:end);
        tagindx = find(strcmpi(taglist(:),prop)); % Find main field
        if ~isempty(tagindx)    % if main field
            eval(['t.handle(i).' taglist{tagindx} ' = val;']);
        else                    % if subfield
            subsindx = findstr(prop,'.');
            if ~isempty(subsindx)
                subtag = prop(1:subsindx-1); % Get main-field
                subsubs = prop(subsindx+1:end);% Get subfield
                ntag = find(strcmpi(taglist(:),subtag));            
                if ~isempty(ntag),                
                    subtagstruct = getfield(t.handle(i),taglist{ntag}); %get substruct
                    if isstruct(subtagstruct)
                        substags = fieldnames(subtagstruct);            
                        nsubs = find(strcmpi(substags(:),subsubs));
                        if ~isempty(nsubs),
                            eval(['t.handle(i).' taglist{ntag} '.' substags{nsubs} ' = val;']);
                        else
                            error('Property Name is not a valid field of CHANHANDLE object');
                        end;
                    else
                        error('Property Name is not a valid field of CHANHANDLE object');
                    end;
                else
                    error('Property Name is not a valid field of CHANHANDLE object');
                end;
            else
                error('Property Name is not a valid field of CHANHANDLE object');
            end;
        end;
    end;
end;