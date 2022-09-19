function [data, mmap] = loadbj(fname,varargin)
%
% data=loadbj(fname,opt)
%    or
% [data, mmap]=loadbj(fname,'param1',value1,'param2',value2,...)
%
% Parse a Binary JData (BJData v1 Draft-2, defined in https://github.com/NeuroJSON/bjdata) 
% file or memory buffer and convert into a MATLAB data structure
%
% By default, this function parses BJD-compliant output. The BJD
% specification is largely similar to UBJSON, with additional data types
% including uint16(u), uint32(m), uint64(M) and half-precision float (h).
% Starting from BJD Draft-2 (JSONLab 3.0 beta or later), all integer and
% floating-point numbers are parsed in Little-Endian as opposed to
% Big-Endian form as in BJD Draft-1/UBJSON Draft-12 (JSONLab 2.0 or older)
%
% authors:Qianqian Fang (q.fang <at> neu.edu)
% initially created on 2013/08/01
%
% input:
%      fname: input file name, if fname contains "{}" or "[]", fname
%             will be interpreted as a BJData/UBJSON string
%      opt: a struct to store parsing options, opt can be replaced by 
%           a list of ('param',value) pairs - the param string is equivallent
%           to a field in opt. opt can have the following 
%           fields (first in [.|.] is the default)
%
%           SimplifyCell [1|0]: if set to 1, loadbj will call cell2mat
%                         for each element of the JSON data, and group 
%                         arrays based on the cell2mat rules.
%           Endian ['L'|'B']: specify the endianness of the numbers
%                         in the BJData/UBJSON input data. Default: 'L'.
%
%                         Starting from JSONLab 2.9, BJData by default uses
%                         [L] Little-Endian for both integers and floating
%                         point numbers. This is a major departure from the
%                         UBJSON specification, where 'B' - Big-Endian -
%                         format is used for integer fields. UBJSON does
%                         not specifically define Endianness for
%                         floating-point numbers, resulting in mixed
%                         implementations. JSONLab 2.0-2.1 used 'B' for
%                         integers and floating-points; JSONLab 1.x uses
%                         'B' for integers and native-endianness for
%                         floating-point numbers.
%           NameIsString [0|1]: for UBJSON Specification Draft 8 or 
%                         earlier versions (JSONLab 1.0 final or earlier), 
%                         the "name" tag is treated as a string. To load 
%                         these UBJSON data, you need to manually set this 
%                         flag to 1.
%           UseMap [0|1]: if set to 1, loadbj uses a containers.Map to 
%                         store map objects; otherwise use a struct object
%           ObjectID [0|interger or list]: if set to a positive number, 
%                         it returns the specified JSON object by index 
%                         in a multi-JSON document; if set to a vector,
%                         it returns a list of specified objects.
%           FormatVersion [2|float]: set the JSONLab format version; since
%                         v2.0, JSONLab uses JData specification Draft 1
%                         for output format, it is incompatible with all
%                         previous releases; if old output is desired,
%                         please set FormatVersion to 1.9 or earlier.
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
%      mmap: (optional) a cell array in the form of
%           {{jsonpath1,[start,length]}, {jsonpath2,[start,length]}, ...}
%           where jsonpath_i is a string in the form of JSONPath, and
%           start is an integer referring to the offset from the begining
%           of the stream, and length is the JSON object string length.
%           For more details, please see the help section of loadjson.m
%
% examples:
%      obj=struct('string','value','array',[1 2 3]);
%      ubjdata=savebj('obj',obj);
%      dat=loadbj(ubjdata)
%      dat=loadbj(['examples' filesep 'example1.bjd'])
%      dat=loadbj(['examples' filesep 'example1.bjd'],'SimplifyCell',0)
%
% license:
%     BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details 
%
% -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%

    opt=varargin2struct(varargin{:});

    if(length(fname)<4096 && exist(fname,'file'))
       fid = fopen(fname,'rb');
       string = fread(fid,jsonopt('MaxBuffer',inf,opt),'uint8=>char')';
       fclose(fid);
    elseif(regexp(fname, '^\s*[\[\{SCHiUIulmLMhdDTFZN]'))
       string=fname;
    else
       error_pos('input file does not exist or buffer is invalid');
    end

    pos = 1;
    inputlen = length(string);
    inputstr = string;

    opt.simplifycell=jsonopt('SimplifyCell',1,opt);
    opt.simplifycellarray=jsonopt('SimplifyCellArray',0,opt);
    opt.usemap=jsonopt('UseMap',0,opt);
    opt.nameisstring=jsonopt('NameIsString',0,opt);
    mmaponly=jsonopt('MmapOnly',0,opt);

    [os,maxelem,systemendian]=computer;
    opt.flipendian_=(systemendian ~= upper(jsonopt('Endian','L',opt)));

    objid=jsonopt('ObjectID',0,opt);
    maxobjid=max(objid);
    if(maxobjid==0)
        maxobjid=inf;
    end

    opt.jsonpath_='$';
    if(nargout>1 || mmaponly)
        mmap={};
    end
    jsoncount=1;
    while pos <= inputlen
        [cc, pos]=next_char(inputstr, pos);
        switch(cc)
            case '{'
                if(nargout>1 || mmaponly)
                    mmap{end+1}={opt.jsonpath_,pos};
                    [data{jsoncount}, pos, newmmap] = parse_object(inputstr, pos, opt);
                    mmap{end}{2}=[mmap{end}{2},pos-mmap{end}{2}];
                    mmap=[mmap(:);newmmap(:)];
                else
                    [data{jsoncount}, pos] = parse_object(inputstr, pos, opt);
                end
            case '['
                if(nargout>1 || mmaponly)
                    mmap{end+1}={opt.jsonpath_,pos};
                    [data{jsoncount}, pos, newmmap] = parse_array(inputstr, pos, opt);
                    mmap{end}{2}=[mmap{end}{2},pos-mmap{end}{2}];
                    mmap=[mmap(:);newmmap(:)];
                else
                    [data{jsoncount}, pos] = parse_array(inputstr, pos, opt);
                end
            case {'S','C','H','i','U','I','u','l','m','L','M','h','d','D','T','F','Z','N'}
                [data{jsoncount}, pos] = parse_value(inputstr, pos, [], opt);
            otherwise
                error_pos('Root level structure must start with a valid marker "{[SCHiUIulmLMhdDTFZN"', inputstr, pos);
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
    if(nargout>1 || mmaponly)
        mmap=mmap';
        mmap=filterjsonmmap(mmap, jsonopt('MMapExclude',{},opt), 0);
        mmap=filterjsonmmap(mmap, jsonopt('MMapInclude',{},opt), 1);
    end
    if(jsonopt('JDataDecode',1,varargin{:})==1)
        try
            data=jdatadecode(data,'Base64',0,'Recursive',1,varargin{:});
        catch ME
            warning(['Failed to decode embedded JData annotations, '...
                'return raw JSON data\n\njdatadecode error: %s\n%s\nCall stack:\n%s\n'], ...
                ME.identifier, ME.message, savejson('',ME.stack));
        end
    end
    if(mmaponly)
        data=mmap;
    end
end

%%-------------------------------------------------------------------------
%% helper functions
%%-------------------------------------------------------------------------

function [data, adv]=parse_block(inputstr, pos, type,count,varargin)
    if(count>=0 && ~isempty(type) && isempty(strfind('iUIulmLMdDh', type)))
        adv=0;
        switch(type)
            case {'S', 'H', '{', '['}
                data=cell(1,count);
                adv=pos;
                for i=1:count
                    [data{i}, pos] = parse_value(inputstr, pos, type, varargin{:});
                end
                adv=pos-adv;
            case 'C'
                data=inputstr(pos:pos+count);
                adv=count;
            case {'T','F','N'}
                error_pos(sprintf('For security reasons, optimized type %c is disabled at position %%d', type),inputstr, pos);
            otherwise
                error_pos(sprintf('Unsupported optimized type %c at position %%d', type),inputstr, pos);
        end
        return;
    end
    [cid,len]=elem_info(inputstr, pos, type);
    datastr=inputstr(pos:pos+len*count-1);
    newdata=uint8(datastr);
    if(varargin{1}.flipendian_)
        newdata=swapbytes(typecast(newdata,cid));
    end
    data=typecast(newdata,cid);
    adv=double(len*count);
end
%%-------------------------------------------------------------------------

function [object, pos, mmap] = parse_array(inputstr,  pos, varargin) % JSON array is written in row-major order
    if(nargout>2)
        mmap={};
        origpath=varargin{1}.jsonpath_;
    end
    pos=parse_char(inputstr, pos, '[');
    object = cell(0, 1);
    dim=[];
    type='';
    count=-1;
    [cc,pos]=next_char(inputstr,pos);
    if(cc == '$')
        type=inputstr(pos+1);
        pos=pos+2;
    end
    [cc,pos]=next_char(inputstr,pos);
    if(cc == '#')
        pos=pos+1;
        [cc,pos]=next_char(inputstr,pos);
        if(cc=='[')
            if(isfield(varargin{1},'noembedding_') && varargin{1}.noembedding_==1)
                error_pos('ND array size specifier does not support embedding');
            end
            varargin{1}.noembedding_=1;
            [dim, pos]=parse_array(inputstr, pos, varargin{:});
            count=prod(double(dim));
            varargin{1}.noembedding_=0;
        else
            [val,pos]=parse_number(inputstr,pos, varargin{:});
            count=double(val);
        end
    end
    if(~isempty(type))
        if(count>=0)
            [object, adv]=parse_block(inputstr, pos, type,count,varargin{:});
            if(~isempty(dim))
                object=permute(reshape(object,fliplr(dim(:)')),length(dim):-1:1);
            end
            pos=pos+adv;
            return;
        else
            endpos=match_bracket(inputstr,pos);
            [cid,len]=elem_info(inputstr, pos, type);
            count=(endpos-pos)/len;
            [object, adv]=parse_block(inputstr, pos, type,count,varargin{:});
            pos=pos+adv;
            pos=parse_char(inputstr, pos, ']');
            return;
        end
    end
    [cc,pos]=next_char(inputstr,pos);
    if cc ~= ']'
         while 1
            if(nargout>2)
                varargin{1}.jsonpath_=[origpath '.' sprintf('[%d]',length(object))];
                mmap{end+1}={varargin{1}.jsonpath_, pos};
                [val, pos, newmmap] = parse_value(inputstr, pos, [], varargin{:});
                mmap{end}{2}=[mmap{end}{2}, pos-mmap{end}{2}];
                mmap=[mmap(:);newmmap(:)];
            else
                [val, pos] = parse_value(inputstr, pos, [], varargin{:});
            end
            object{end+1} = val;
            [cc,pos]=next_char(inputstr,pos);
            if cc == ']'
                break;
            end
         end
    end
    if(varargin{1}.simplifycell)
      if(iscell(object) && ~isempty(object) && isnumeric(object{1}))
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
    if(count==-1)
        pos=parse_char(inputstr, pos, ']');
    end
end

%%-------------------------------------------------------------------------

function pos=parse_char(inputstr, pos, c)
    if pos > length(inputstr) || inputstr(pos) ~= c
        error_pos(sprintf('Expected %c at position %%d', c),inputstr, pos);
    else
        pos = pos + 1;
    end
end

%%-------------------------------------------------------------------------

function [c, pos] = next_char(inputstr, pos)
    if pos > length(inputstr)
        c = [];
    else
        c = inputstr(pos);
        while(c=='N')
             pos=pos+1;
             c = inputstr(pos);
        end
    end
end

%%-------------------------------------------------------------------------
function [str, pos] = parse_name(inputstr, pos, varargin)
    [val, pos]=parse_number(inputstr,pos,varargin{:});
    bytelen=double(val);
    if(length(inputstr)>=pos+bytelen-1)
        str=inputstr(pos:pos+bytelen-1);
        pos=pos+bytelen;
    else
        error_pos('End of file while expecting end of name', inputstr, pos);
    end
end

%%-------------------------------------------------------------------------

function [str, pos] = parseStr(inputstr, pos, type, varargin)
    if(isempty(type))
        type=inputstr(pos);
        if type ~= 'S' && type ~= 'C' && type ~= 'H'
            error_pos('String starting with S expected at position %d',inputstr, pos);
        else
            pos = pos + 1;
        end
    end

    if(type == 'C')
        str=inputstr(pos);
        pos=pos+1;
        return;
    end
    [val, pos]=parse_number(inputstr,pos,varargin{:});
    bytelen=double(val);
    if(length(inputstr)>=pos+bytelen-1)
        str=inputstr(pos:pos+bytelen-1);
        pos=pos+bytelen;
    else
        error_pos('End of file while expecting end of inputstr',inputstr, pos);
    end
end

%%-------------------------------------------------------------------------

function [num, pos] = parse_number(inputstr, pos, varargin)
    id=strfind('iUIulmLMhdD',inputstr(pos));
    if(isempty(id))
        error_pos('expecting a number at position %d',inputstr, pos);
    end
    type={'int8','uint8','int16','uint16','int32','uint32','int64','uint64','half','single','double'};
    bytelen=[1,1,2,2,4,4,8,8,2,4,8];
    if(~exist('half','builtin'))
        type{9}='uint16';
    end
    datastr=inputstr(pos+1:pos+bytelen(id));
    newdata=uint8(datastr);
    if(varargin{1}.flipendian_)
        newdata=swapbytes(typecast(newdata,type{id}));
    end
    num=typecast(newdata,type{id});
    pos = pos + bytelen(id)+1;
end

%%-------------------------------------------------------------------------

function varargout = parse_value(inputstr, pos, type, varargin)
    if(length(type)==1)
            cc=type;
            varargout{2}=pos;
    else
            [cc,varargout{2}]=next_char(inputstr,pos);
    end
    if(nargout>2)
            varargout{3}={};
    end
    switch(cc)
        case {'S','C','H'}
            [varargout{1:2}] = parseStr(inputstr, varargout{2}, type, varargin{:});
            return;
        case '['
            [varargout{1:nargout}] = parse_array(inputstr, varargout{2}, varargin{:});
            return;
        case '{'
            [varargout{1:nargout}] = parse_object(inputstr, varargout{2}, varargin{:});
            return;
        case {'i','U','I','u','l','m','L','M','h','d','D'}
            [varargout{1:2}] = parse_number(inputstr, varargout{2}, varargin{:});
            return;
        case 'T'
            varargout{1} = true;
            varargout{2} = varargout{2} + 1;
            return;
        case 'F'
            varargout{1} = false;
            varargout{2} = varargout{2} + 1;
            return;
        case {'Z','N'}
            varargout{1} = [];
            varargout{2} = varargout{2} + 1;
            return;
    end
    error_pos('Value expected at position %d', inputstr, varargout{2});
end
%%-------------------------------------------------------------------------

function pos=error_pos(msg, inputstr, pos)
    poShow = max(min([pos-15 pos-1 pos pos+20],length(inputstr)),1);
    if poShow(3) == poShow(2)
        poShow(3:4) = poShow(2)+[0 -1];  % display nothing after
    end
    msg = [sprintf(msg, pos) ': ' ...
    inputstr(poShow(1):poShow(2)) '<error>' inputstr(poShow(3):poShow(4)) ];
    error( ['JSONLAB:BJData:InvalidFormat: ' msg] );
end

%%-------------------------------------------------------------------------
function [object, pos, mmap] = parse_object(inputstr, pos, varargin)
    if(nargout>2)
        mmap={};
        origpath=varargin{1}.jsonpath_;
    end
    pos=parse_char(inputstr,pos,'{');
    usemap=varargin{1}.usemap;
    if(usemap)
	object = containers.Map();
    else
	object = [];
    end
    count=-1;
    [cc, pos]=next_char(inputstr,pos);
    if(cc == '$')
        pos=pos+2;
    end
    [cc, pos]=next_char(inputstr,pos);
    if(cc == '#')
        pos=pos+1;
        [val,pos]=parse_number(inputstr, pos, varargin{:});
        count=double(val);
    end
    [cc, pos]=next_char(inputstr,pos);
    if cc ~= '}'
        num=0;
        while 1
            if(varargin{1}.nameisstring)
                [str, pos] = parseStr(inputstr, pos, [], varargin{:});
            else
                [str, pos] = parse_name(inputstr, pos, varargin{:});
            end
            if isempty(str)
                str='x0x0_'; % empty name is valid in BJData/UBJSON, decodevarname('x0x0_') restores '\0'
            end
            if(nargout>2)
                varargin{1}.jsonpath_=[origpath,'.',str];
                mmap{end+1}={varargin{1}.jsonpath_,pos};
                [val, pos,newmmap] = parse_value(inputstr, pos, [], varargin{:});
                mmap{end}{2}=[mmap{end}{2}, pos-mmap{end}{2}];
                mmap=[mmap(:);newmmap(:)];
            else
                [val, pos] = parse_value(inputstr, pos, [], varargin{:});
            end
            num=num+1;
            if(usemap)
                object(str)=val;
            else
                object.(encodevarname(str,varargin{:}))=val;
            end
            [cc, pos]=next_char(inputstr,pos);
            if cc == '}' || (count>=0 && num>=count)
                break;
            end
        end
    end
    if(count==-1)
        pos=parse_char(inputstr, pos, '}');
    end
end

%%-------------------------------------------------------------------------
function [cid,len]=elem_info(inputstr, pos, type)
    id=strfind('iUIulmLMhdD',type);
    type={'int8','uint8','int16','uint16','int32','uint32','int64','uint64','half','single','double'};
    bytelen=[1,1,2,2,4,4,8,8,2,4,8];
    if(~exist('half','builtin'))
        type{9}='uint16';
    end
    if(id>0)
        cid=type{id};
        len=bytelen(id);
    else
        error_pos('unsupported type at position %d',inputstr, pos);
    end
end
