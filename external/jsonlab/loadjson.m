function [data, mmap] = loadjson(fname,varargin)
%
% data=loadjson(fname,opt)
%    or
% [data, mmap]=loadjson(fname,'param1',value1,'param2',value2,...)
%
% parse a JSON (JavaScript Object Notation) file or string and return a
% matlab data structure with optional memory-map (mmap) table
%
% authors:Qianqian Fang (q.fang <at> neu.edu)
% created on 2011/09/09, including previous works from 
%
%         Nedialko Krouchev: http://www.mathworks.com/matlabcentral/fileexchange/25713
%            created on 2009/11/02
%         François Glineur: http://www.mathworks.com/matlabcentral/fileexchange/23393
%            created on  2009/03/22
%         Joel Feenstra:
%         http://www.mathworks.com/matlabcentral/fileexchange/20565
%            created on 2008/07/03
%
% input:
%      fname: input file name; if fname contains "{}" or "[]", fname
%             will be interpreted as a JSON string
%      opt: (optional) a struct to store parsing options, opt can be replaced by 
%           a list of ('param',value) pairs - the param string is equivallent
%           to a field in opt. opt can have the following 
%           fields (first in [.|.] is the default)
%
%           SimplifyCell [1|0]: if set to 1, loadjson will call cell2mat
%                         for each element of the JSON data, and group 
%                         arrays based on the cell2mat rules.
%           FastArrayParser [1|0 or integer]: if set to 1, use a
%                         speed-optimized array parser when loading an 
%                         array object. The fast array parser may 
%                         collapse block arrays into a single large
%                         array similar to rules defined in cell2mat; 0 to 
%                         use a legacy parser; if set to a larger-than-1
%                         value, this option will specify the minimum
%                         dimension to enable the fast array parser. For
%                         example, if the input is a 3D array, setting
%                         FastArrayParser to 1 will return a 3D array;
%                         setting to 2 will return a cell array of 2D
%                         arrays; setting to 3 will return to a 2D cell
%                         array of 1D vectors; setting to 4 will return a
%                         3D cell array.
%           UseMap [0|1]: if set to 1, loadjson uses a containers.Map to 
%                         store map objects; otherwise use a struct object
%           ShowProgress [0|1]: if set to 1, loadjson displays a progress bar.
%           ParseStringArray [0|1]: if set to 0, loadjson converts "string arrays" 
%                         (introduced in MATLAB R2016b) to char arrays; if set to 1,
%                         loadjson skips this conversion.
%           FormatVersion [3|float]: set the JSONLab format version; since
%                         v2.0, JSONLab uses JData specification Draft 1
%                         for output format, it is incompatible with all
%                         previous releases; if old output is desired,
%                         please set FormatVersion to 1.9 or earlier.
%           Encoding ['']: json file encoding. Support all encodings of
%                         fopen() function
%           ObjectID [0|interger or list]: if set to a positive number, 
%                         it returns the specified JSON object by index 
%                         in a multi-JSON document; if set to a vector,
%                         it returns a list of specified objects.
%           JDataDecode [1|0]: if set to 1, call jdatadecode to decode
%                         JData structures defined in the JData
%                         Specification.
%           BuiltinJSON [0|1]: if set to 1, this function attempts to call
%                         jsondecode, if presents (MATLAB R2016b or Octave
%                         6) first. If jsondecode does not exist or failed, 
%                         this function falls back to the jsonlab parser
%           MmapOnly [0|1]: if set to 1, this function only returns mmap
%           MMapInclude 'str1' or  {'str1','str2',..}: if provided, the
%                         returned mmap will be filtered by only keeping
%                         entries containing any one of the string patterns
%                         provided in a cell
%           MMapExclude 'str1' or  {'str1','str2',..}: if provided, the
%                         returned mmap will be filtered by removing
%                         entries containing any one of the string patterns
%                         provided in a cell
%
% output:
%      dat: a cell array, where {...} blocks are converted into cell arrays,
%           and [...] are converted to arrays
%      mmap: (optional) a cell array as memory-mapping table in the form of
%             {{jsonpath1,[start,length,<whitespace>]},
%              {jsonpath2,[start,length,<whitespace>]}, ...}
%           where jsonpath_i is a string in the JSONPath [1,2] format, and
%           "start" is an integer referring to the offset from the begining
%           of the stream, and "length" is the JSON object string length.
%           An optional 3rd integer "whitespace" may appear to record the
%           preceding whitespace length in case expansion of the data
%           record is needed when using the mmap.
%
%           Memory-mapping table (mmap) is useful when fast reading/writing
%           specific data records inside a large JSON file without needing
%           to load/parse/overwrite the entire file.
%
%           The JSONPath keys used in mmap is largely compatible to the
%           upstream specification defined in [1], with a slight extension
%           to handle contatenated JSON files.
%
%           In the mmap jsonpath key, a '$' denotes the root object, a '.'
%           denotes a child of the preceding element; '.key' points to the
%           value segment of the child named "key" of the preceding
%           object; '.[i]' denotes the (i+1)th member of the preceding
%           element, which must be an array. For example, a key
%
%           $.obj1.obj2.[0].obj3
%
%           defines the memory-map of the "value" section in the below
%           hierarchy:
%             {
%                "obj1":{
%            	     "obj2":[
%                       {"obj3":value},
%                       ...
%                    ],
%                    ...
%                 }
%             }
%           Please note that "value" can be any valid JSON value, including
%           an array, an object, a string or numerical value.
%
%           To handle concatenated JSON objects (including ndjson,
%           http://ndjson.org/), such as
%
%             {"root1": {"obj1": ...}}
%             ["root2", value1, value2, {"obj2": ...}]
%             {"root3": ...}
%
%           we use '$' or '$0' for the first root-object, and '$1' refers
%           to the 2nd root object (["root2",...]) and '$2' referrs to the
%           3rd root object, and so on. Please note that this syntax is an
%           extension from the JSONPath documentation [1,2]
%
%           [1] https://goessner.net/articles/JsonPath/
%           [2] http://jsonpath.herokuapp.com/
%
% examples:
%      dat=loadjson('{"obj":{"string":"value","array":[1,2,3]}}')
%      dat=loadjson(['examples' filesep 'example1.json'])
%      [dat, mmap]=loadjson(['examples' filesep 'example1.json'],'SimplifyCell',0)
%
% license:
%     BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details 
%
% -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%

    opt=varargin2struct(varargin{:});
    
    if(regexp(fname,'^\s*(?:\[.*\])|(?:\{.*\})\s*$','once'))
       string=fname;
    elseif(exist(fname,'file'))
       try
           encoding = jsonopt('Encoding','',opt);
           if(isempty(encoding))
               string = fileread(fname);
           else
               fid = fopen(fname,'r','n',encoding);
               string = fread(fid,'*char')';
               fclose(fid);
           end
       catch
           try
               string = urlread(fname);
           catch
               string = urlread(['file://',fullfile(pwd,fname)]);
           end
       end
    elseif(regexpi(fname,'^\s*(http|https|ftp|file)://'))
       string = urlread(fname);
    else
       error_pos('input file does not exist');
    end

    if(jsonopt('BuiltinJSON',0,opt) && exist('jsondecode','builtin'))
        try
            newstring=regexprep(string,'[\r\n]','');
            newdata=jsondecode(newstring);
            newdata=jdatadecode(newdata,'Base64',1,'Recursive',1,varargin{:});
            data=newdata;
            return;
        catch
            warning('built-in jsondecode function failed to parse the file, fallback to loadjson');
        end
    end

    pos = 1; inputlen = length(string); inputstr = string;
    arraytokenidx=find(inputstr=='[' | inputstr==']');
    arraytoken=inputstr(arraytokenidx);

    % String delimiters and escape chars identified to improve speed:
    esc = find(inputstr=='"' | inputstr=='\' ); % comparable to: regexp(inputstr, '["\\]');
    index_esc = 1;

    opt.arraytoken_=arraytoken;
    opt.arraytokenidx_=arraytokenidx;
    opt.simplifycell=jsonopt('SimplifyCell',1,opt);
    opt.simplifycellarray=jsonopt('SimplifyCellArray',opt.simplifycell,opt);
    opt.formatversion=jsonopt('FormatVersion',3,opt);
    opt.fastarrayparser=jsonopt('FastArrayParser',1,opt);
    opt.parsestringarray=jsonopt('ParseStringArray',0,opt);
    opt.usemap=jsonopt('UseMap',0,opt);
    opt.arraydepth_=1;
    opt.mmaponly=jsonopt('MmapOnly',0,opt);

    if(jsonopt('ShowProgress',0,opt)==1)
        opt.progressbar_=waitbar(0,'loading ...');
    end

    objid=jsonopt('ObjectID',0,opt);
    maxobjid=max(objid);
    if(maxobjid==0)
        maxobjid=inf;
    end
    opt.jsonpath_='$';
    if(nargout>1 || opt.mmaponly)
        mmap={};
    end
    jsoncount=1;
    while pos <= inputlen
        [cc,pos,w1]=next_char(inputstr, pos);
        switch(cc)
            case '{'
                if(nargout>1 || opt.mmaponly)
                    mmap{end+1}={opt.jsonpath_,[pos, 0, w1]};
                    [data{jsoncount},pos,index_esc,newmmap] = parse_object(inputstr, pos, esc, index_esc,opt);
                    mmap{end}{2}(2)=pos-mmap{end}{2}(1);
                    mmap=[mmap(:);newmmap(:)];
                else
                    [data{jsoncount},pos,index_esc] = parse_object(inputstr, pos, esc, index_esc,opt);
                end
            case '['
                if(nargout>1 || opt.mmaponly)
                    mmap{end+1}={opt.jsonpath_,[pos,0,w1]};
                    [data{jsoncount},pos,index_esc,newmmap] = parse_array(inputstr, pos, esc, index_esc,opt);
                    mmap{end}{2}(2)=pos-mmap{end}{2}(1);
                    mmap=[mmap(:);newmmap(:)];
                else
                    [data{jsoncount},pos,index_esc] = parse_array(inputstr, pos, esc, index_esc,opt);
                end
            otherwise
                pos=error_pos('Outer level structure must be an object or an array',inputstr,pos);
        end
	    if(jsoncount>=maxobjid)
	        break;
	    end
        opt.jsonpath_=sprintf('$%d',jsoncount);
        jsoncount=jsoncount+1;
    end % while

    if(length(objid)>1 || min(objid)>1)
        data=data(objid(objid<=length(data)));
    end

    jsoncount=length(data);
    if(jsoncount==1 && iscell(data))
        data=data{1};
    end
    if(nargout>1 || opt.mmaponly)
        mmap=mmap';
        mmap=filterjsonmmap(mmap, jsonopt('MMapExclude',{},opt), 0);
        mmap=filterjsonmmap(mmap, jsonopt('MMapInclude',{},opt), 1);
        mmap=cellfun(@(x) {x{1},x{2}(1:(2+int8(length(x{2})>=3 && (x{2}(3)>0))))}, mmap, 'UniformOutput', false);
    end
    if(jsonopt('JDataDecode',1,varargin{:})==1)
        try
            data=jdatadecode(data,'Base64',1,'Recursive',1,varargin{:});
        catch ME
            warning(['Failed to decode embedded JData annotations, '...
                'return raw JSON data\n\njdatadecode error: %s\n%s\nCall stack:\n%s\n'], ...
                ME.identifier, ME.message, savejson('',ME.stack));
        end
    end
    if(opt.mmaponly)
        data=mmap;
    end
    if(isfield(opt,'progressbar_'))
        close(opt.progressbar_);
    end
end

%%-------------------------------------------------------------------------
%% helper functions
%%-------------------------------------------------------------------------

function [object, pos,index_esc, mmap] = parse_array(inputstr, pos, esc, index_esc, varargin) % JSON array is written in row-major order
    if(nargout>3)
        mmap={};
        origpath=varargin{1}.jsonpath_;
    end
    pos=parse_char(inputstr, pos, '[');
    object = cell(0, 1);
    arraydepth=varargin{1}.arraydepth_;
    pbar=-1;
    if(isfield(varargin{1},'progressbar_'))
        pbar=varargin{1}.progressbar_;
    end
    format=varargin{1}.formatversion;
    [cc,pos]=next_char(inputstr,pos);
    endpos=[];
    
    if cc ~= ']'
        try
            if((varargin{1}.fastarrayparser)>=1 && arraydepth>=varargin{1}.fastarrayparser)
                [endpos, maxlevel]=fast_match_bracket(varargin{1}.arraytoken_,varargin{1}.arraytokenidx_,pos);
                if(~isempty(endpos))
                    arraystr=['[' inputstr(pos:endpos)];
                    arraystr=sscanf_prep(arraystr);
                    if(isempty(find(arraystr=='"', 1)))
                        % handle 1D array first
                        if(maxlevel==1)
                            astr=arraystr(2:end-1);
                            astr(astr==' ')='';
                            [obj, count, errmsg, nextidx]=sscanf(astr,'%f,',[1,inf]);
                            if(nextidx>=length(astr)-1)
                                    object=obj;
                                    pos=endpos;
                                    pos=parse_char(inputstr, pos, ']');
                                    return;
                            end
                        end

                        % for N-D packed array in a nested array construct, 
                        if(maxlevel>=2 && ~isempty(regexp(arraystr(2:end),'^\s*\[','once')))
                            [dims,isndarray]=nestbracket2dim(arraystr);
                            rowstart=find(arraystr(2:end)=='[',1)+1;
                            if(rowstart && isndarray)
                                [obj, nextidx]=parsendarray(arraystr,dims);
                                if(nextidx>=length(arraystr)-1)
                                    object=obj;
                                    if(format>1.9)
                                        object=permute(object,ndims(object):-1:1);
                                    end
                                    pos=endpos;
                                    pos=parse_char(inputstr, pos, ']');
                                    if(pbar>0)
                                        waitbar(pos/length(inputstr),pbar,'loading ...');
                                    end
                                    return;
                                end
                            end
                        end
                    end
                end
            end
            if(isempty(regexp(arraystr,':','once')) && isempty(regexp(arraystr,'\(','once')))
                arraystr=regexprep(arraystr,'\[','{');
                arraystr=regexprep(arraystr,'\]','}');
                if(varargin{1}.parsestringarray==0)
                    arraystr=regexprep(arraystr,'\"','''');
                end
                object=eval(arraystr);
                if(iscell(object))
                    object=cellfun(@unescapejsonstring,object,'UniformOutput',false);
                end
                pos=endpos;
            end
        catch
        end
        if(isempty(endpos) || pos~=endpos)
            w2=0;
            while 1
                varargin{1}.arraydepth_=arraydepth+1;
                if(nargout>3)
                    varargin{1}.jsonpath_=[origpath '.' sprintf('[%d]',length(object))];
                    mmap{end+1}={varargin{1}.jsonpath_, [pos, 0, w2]};
                    [val, pos, index_esc, newmmap] = parse_value(inputstr, pos, esc, index_esc,varargin{:});
                    mmap{end}{2}(2)=pos-mmap{end}{2}(1);
                    mmap=[mmap(:);newmmap(:)];
                else
                    [val, pos,index_esc] = parse_value(inputstr, pos, esc, index_esc,varargin{:});
                end
                object{end+1} = val;
                [cc,pos]=next_char(inputstr,pos);
                if cc == ']'
                    break;
                end
                [pos, w1, w2]=parse_char(inputstr, pos, ',');
            end
        end
    end

    if(varargin{1}.simplifycell)
      if(iscell(object) && ~isempty(object) && (isnumeric(object{1}) || isstruct(object{1})) )
          if(all(cellfun(@(e) isequal(size(object{1}), size(e)) , object(2:end))))
              try
                  oldobj=object;
                  if(iscell(object) && length(object)>1 && ndims(object{1})>=2)
                      catdim=size(object{1});
                      catdim=ndims(object{1})-(catdim(end)==1)+1;
                      object=cat(catdim,object{:});
                      object=permute(object,ndims(object):-1:1);
                  else
                      object=cell2mat(object')';
                  end
                  if(iscell(oldobj) && isstruct(object) && numel(object)>1 && varargin{1}.simplifycellarray==0)
                      object=oldobj;
                  end
              catch
              end
          end
      end
      if(~iscell(object) && size(object,1)>1 && ndims(object)==2)
            object=object';
      end
    end
    pos=parse_char(inputstr, pos, ']');

    if(pbar>0)
        waitbar(pos/length(inputstr),pbar,'loading ...');
    end
end
%%-------------------------------------------------------------------------

function [pos, w1, w2]=parse_char(inputstr, pos, c)
    w1=pos;
    w2=0;
    pos=skip_whitespace(pos, inputstr);
    w1=pos-w1;
    if pos > length(inputstr) || inputstr(pos) ~= c
        pos=error_pos(sprintf('Expected %c at position %%d', c),inputstr,pos);
    else
        pos = pos + 1;
        w2=pos;
        pos=skip_whitespace(pos, inputstr);
        w2=pos-w2;
    end
end
%%-------------------------------------------------------------------------

function [c, pos, w1] = next_char(inputstr, pos)
    w1=pos;
    pos=skip_whitespace(pos, inputstr);
    w1=pos-w1;
    if pos > length(inputstr)
        c = [];
    else
        c = inputstr(pos);
    end
end

%%-------------------------------------------------------------------------
function [str, pos,index_esc, mmap] = parseStr(inputstr, pos, esc, index_esc, varargin)
    if(nargout>3)
        mmap={};
    end
    if inputstr(pos) ~= '"'
        pos=error_pos('String starting with " expected at position %d',inputstr,pos);
    else
        pos = pos + 1;
    end
    str = '';
    while pos <= length(inputstr)
        while index_esc <= length(esc) && esc(index_esc) < pos
            index_esc = index_esc + 1;
        end
        if index_esc > length(esc)
            str = [str inputstr(pos:end)];
            pos = length(inputstr) + 1;
            break;
        else
            str = [str inputstr(pos:esc(index_esc)-1)];
            pos = esc(index_esc);
        end
        nstr = length(str);
        switch inputstr(pos)
            case '"'
                pos = pos + 1;
                if(~isempty(str))
                    if(strcmp(str,'_Inf_'))
                        str=Inf;
                    elseif(strcmp(str,'-_Inf_'))
                        str=-Inf;
                    elseif(strcmp(str,'_NaN_'))
                        str=NaN;
                    end
                end
                return;
            case '\'
                if pos+1 > length(inputstr)
                    pos=error_pos('End of file reached right after escape character',inputstr,pos);
                end
                pos = pos + 1;
                switch inputstr(pos)
                    case {'"' '\' '/'}
                        str(nstr+1) = inputstr(pos);
                        pos = pos + 1;
                    case {'b' 'f' 'n' 'r' 't'}
                        str(nstr+1) = sprintf(['\' inputstr(pos)]);
                        pos = pos + 1;
                    case 'u'
                        if pos+4 > length(inputstr)
                            pos=error_pos('End of file reached in escaped unicode character',inputstr,pos);
                        end
                        str(nstr+(1:6)) = inputstr(pos-1:pos+4);
                        pos = pos + 5;
                end
            otherwise % should never happen
                str(nstr+1) = inputstr(pos);
                keyboard;
                pos = pos + 1;
        end
    end
    str=unescapejsonstring(str);
    pos=error_pos('End of file while expecting end of inputstr',inputstr,pos);
end
%%-------------------------------------------------------------------------

function [num, pos] = parse_number(inputstr, pos, varargin)
    currstr=inputstr(pos:min(pos+30,end));
    [num, one, err, delta] = sscanf(currstr, '%f', 1);
    if ~isempty(err)
        pos=error_pos('Error reading number at position %d',inputstr,pos);
    end
    pos = pos + delta-1;
end
%%-------------------------------------------------------------------------

function varargout = parse_value(inputstr, pos, esc, index_esc, varargin)
    len=length(inputstr);
    if(isfield(varargin{1},'progressbar_'))
        waitbar(pos/len,varargin{1}.progressbar_,'loading ...');
    end
    varargout{3}=index_esc;
    if(nargout>3)
            varargout{4}={};
    end
    switch(inputstr(pos))
        case '"'
            [varargout{1:nargout}] = parseStr(inputstr, pos, esc, index_esc,varargin{:});
            return;
        case '['
            [varargout{1:nargout}] = parse_array(inputstr, pos, esc, index_esc, varargin{:});
            return;
        case '{'
            [varargout{1:nargout}] = parse_object(inputstr, pos, esc, index_esc, varargin{:});
            return;
        case {'-','0','1','2','3','4','5','6','7','8','9'}
            [varargout{1:2}] = parse_number(inputstr, pos, varargin{:});
            return;
        case 't'
            if pos+3 <= len && strcmpi(inputstr(pos:pos+3), 'true')
                varargout{1} = true;
                varargout{2} = pos + 4;
                return;
            end
        case 'f'
            if pos+4 <= len && strcmpi(inputstr(pos:pos+4), 'false')
                varargout{1} = false;
                varargout{2} = pos + 5;
                return;
            end
        case 'n'
            if pos+3 <= len && strcmpi(inputstr(pos:pos+3), 'null')
                varargout{1} = [];
                varargout{2} = pos + 4;
                return;
            end
    end
    varargout{2}=error_pos('Value expected at position %d',inputstr,pos);
end

%%-------------------------------------------------------------------------
function [object, pos, index_esc, mmap] = parse_object(inputstr, pos, esc, index_esc, varargin)
    if(nargout>3)
        mmap={};
        origpath=varargin{1}.jsonpath_;
    end
    pos=parse_char(inputstr, pos, '{');
    usemap=varargin{1}.usemap;
    if(usemap)
	object = containers.Map();
    else
	object = [];
    end
    [cc,pos]=next_char(inputstr,pos);
    if cc ~= '}'
        while 1
            [str, pos, index_esc] = parseStr(inputstr, pos, esc, index_esc, varargin{:});
            if isempty(str) && ~usemap
                str='x0x0_'; % empty name is valid in JSON, decodevarname('x0x0_') restores '\0'
            end
            [pos, w1, w2]=parse_char(inputstr, pos, ':');
            if(nargout>3)
                varargin{1}.jsonpath_=[origpath,'.',str];
                mmap{end+1}={varargin{1}.jsonpath_,[pos,0,w2]};
                [val, pos,index_esc, newmmap] = parse_value(inputstr, pos, esc, index_esc, varargin{:});
                mmap{end}{2}(2)=pos-mmap{end}{2}(1);
                mmap=[mmap(:);newmmap(:)];
            else
                [val, pos,index_esc] = parse_value(inputstr, pos, esc, index_esc, varargin{:});
            end
            if(usemap)
		object(str)=val;
	    else
		object.(encodevarname(str,varargin{:}))=val;
	    end
            [cc,pos]=next_char(inputstr,pos);
            if cc == '}'
                break;
            end
            pos=parse_char(inputstr, pos, ',');
        end
    end
    pos=parse_char(inputstr, pos, '}');
end

%%-------------------------------------------------------------------------

function pos=error_pos(msg, inputstr, pos)
    poShow = max(min([pos-15 pos-1 pos pos+20],length(inputstr)),1);
    if poShow(3) == poShow(2)
        poShow(3:4) = poShow(2)+[0 -1];  % display nothing after
    end
    msg = [sprintf(msg, pos) ': ' ...
    inputstr(poShow(1):poShow(2)) '<error>' inputstr(poShow(3):poShow(4)) ];
    error( ['JSONLAB:JSON:InvalidFormat: ' msg] );
end

%%-------------------------------------------------------------------------

function newpos=skip_whitespace(pos, inputstr)
    newpos=pos;
    while newpos <= length(inputstr) && isspace(inputstr(newpos))
        newpos = newpos + 1;
    end
end

%%-------------------------------------------------------------------------
function newstr=unescapejsonstring(str)
    newstr=str;
    if(iscell(str))
        try
            newstr=cell2mat(cellfun(@(x) cell2mat(x),str(:),'un',0));
        catch
        end
    end
    if(~ischar(str))
        return;
    end
    escapechars={'\\','\"','\/','\a','\b','\f','\n','\r','\t','\v'};
    for i=1:length(escapechars)
        newstr=regexprep(newstr,regexprep(escapechars{i},'\\','\\\\'), escapechars{i});
    end
    newstr=regexprep(newstr,'\\u([0-9A-Fa-f]{4})', '${char(base2dec($1,16))}');
end

%%-------------------------------------------------------------------------
function arraystr=sscanf_prep(str)
    arraystr=str;
    if(regexp(str,'"','once'))
        arraystr=regexprep(arraystr,'"_NaN_"','NaN');
        arraystr=regexprep(arraystr,'"([-+]*)_Inf_"','$1Inf');
    end
    arraystr(arraystr==sprintf('\n'))=' ';
    arraystr(arraystr==sprintf('\r'))=' ';
end

%%-------------------------------------------------------------------------
function [obj, nextidx]=parsendarray(arraystr, dims)
    astr=arraystr;
    astr(astr=='[')=' ';
    astr(astr==']')=' ';
    astr(astr==',')=' ';
    [obj, count, errmsg, nextidx]=sscanf(astr,'%f',inf);
    if(nextidx>=length(astr)-1)
            obj=reshape(obj,dims);
            nextidx=length(arraystr)+1;
    end
end
