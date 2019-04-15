function t = set(t,varargin)
% T = SET(T.PROP,VAL); Set object properties
% Usage:
% t = set(t,property1,value1,property2,value2,...)
%	Sets values of properties 'property1', 'property2', etc. to values 
%	'value1', value2', ... respectively

if nargin<3
    error('Sorry, Insufficient number of arguments');
elseif rem(nargin-1,2)
    error('Sorry, Property and value must be given in pairs');
else
    if ~isa(t,'acqparam')
        error('First input must belong to class acqparam');
    end;
    inargs = varargin;              % Input arguments
    taglist = fieldnames(struct(t)); % List of fields
    while length(inargs)>=2
        prop = inargs{1};           % One property at a time
        val  = inargs{2};
        inargs = inargs(3:end);
        tagindx = find(strcmpi(taglist(:),prop)); % Find which property?
        if ~isempty(tagindx)        % If regular field
            eval(['t.' taglist{tagindx} ' = val;']);
        else                        % If SampleInfo field
            sampinfotaglist = fieldnames(t.SampleInfo);
            subindx = find(strcmpi(sampinfotaglist(:),prop));
            if ~isempty(subindx)
                eval(['t.SampleInfo.' sampinfotaglist{subindx} '=val;']);
            else
                error(['Couldnt find requested field - ''',prop, ''' to update']);
            end;
        end;
    end;
end;
