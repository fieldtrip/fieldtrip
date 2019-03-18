function t = set(t,varargin)
% T = SET(T,PROP,VAL); Set object (T) property (PROP) to specified value (VAL)
% Usage:
% t = set(t,property1,value1,property2,value2,...)
%	Sets values of properties 'property1', 'property2', etc. to values 
%	'value1', value2', ... respectively

switch nargin
case {0,1,2}
    error('Sorry, Atleast three inputs are required');
otherwise
    if ~isa(t,'sqdhandle')
        error('First input must belong to class sqdhandle');
    end;
    if rem(nargin-1,2)
        error('Property and values must be given as pairs');
    end;
    inargs = varargin;              % List o input arguments
    taglist = fieldnames(struct(t)); % Fields of objects
    while length(inargs)>=2
        prop = inargs{1};           % One property at a time
        val  = inargs{2};
        inargs = inargs(3:end);
        tagindx = find(strcmpi(taglist(:),prop));
        if ~isempty(tagindx)        % If main field of object
            eval(['t.' taglist{tagindx} ' = val;']);
            if strcmpi(taglist{tagindx},'FileName')  % If new sqd-file
                t = readsqdinfo(t,get(t,'FileName')); % update all fields
            end;
        else
            patinfolist = fieldnames(t.PatientInfo); % fields of patient info
            subindx = find(strcmpi(patinfolist(:),prop));
            if ~subindx             % if field of patient info
                eval(['t.PatientInfo.' patinfolist{subindx} ' = val;']);
            else                    % if acqparam
                t.acqparam = set(t.acqparam,prop,val);
            end;            
        end;
    end;
end;