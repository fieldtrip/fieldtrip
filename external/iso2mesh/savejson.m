function json=savejson(rootname,obj,varargin)
%
% json=savejson(rootname,obj,opt)
%    or
% json=savejson(rootname,obj,'param1',value1,'param2',value2,...)
%
% convert a MATLAB object (cell, struct or array) into a JSON (JavaScript
% Object Notation) string
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%            created on 2011/09/09
%
% $Id$
%
% input:
%      rootname: name of the root-object, if set to '', will use variable name
%      obj: a MATLAB object (array, cell, cell array, struct, struct array)
%      opt: a struct for additional options, use [] if all use default
%        opt can have the following fields (first in [.|.] is the default)
%        opt.FloatFormat ['%.10g'|string]: format to show each numeric element
%                         of a 1D/2D array;
%        opt.ArrayIndent [1|0]: if 1, output explicit data array with
%                         precedent indentation; if 0, no indentation
%        opt.ArrayToStruct[0|1]: when set to 0, savejson outputs 1D/2D
%                         array in JSON array format; if sets to 1, an
%                         array will be shown as a struct with fields
%                         "_ArrayType_", "_ArraySize_" and "_ArrayData_"; for
%                         sparse arrays, the non-zero elements will be
%                         saved to _ArrayData_ field in triplet-format i.e.
%                         (ix,iy,val) and "_ArrayIsSparse_" will be added
%                         with a value of 1; for a complex array, the 
%                         _ArrayData_ array will include two columns 
%                         (4 for sparse) to record the real and imaginary 
%                         parts, and also "_ArrayIsComplex_":1 is added. 
%        opt.ParseLogical [0|1]: if this is set to 1, logical array elem
%                         will use true/false rather than 1/0.
%        opt.NoRowBracket [1|0]: if this is set to 1, arrays with a single
%                         numerical element will be shown without a square
%                         bracket, unless it is the root object; if 0, square
%                         brackets are forced for any numerical arrays.
%        opt.ForceRootName [0|1]: when set to 1 and rootname is empty, savejson
%                         will use the name of the passed obj variable as the 
%                         root object name; if obj is an expression and 
%                         does not have a name, 'root' will be used; if this 
%                         is set to 0 and rootname is empty, the root level 
%                         will be merged down to the lower level.
%        opt can be replaced by a list of ('param',value) pairs. The param 
%        string is equivallent to a field in opt.
% output:
%      json: a string in the JSON format (see http://json.org)
%
% examples:
%      a=struct('node',[1  9  10; 2 1 1.2], 'elem',[9 1;1 2;2 3],...
%           'face',[9 01 2; 1 2 3; NaN,Inf,-Inf], 'author','FangQ');
%      savejson('mesh',a)
%      savejson('',a,'ArrayIndent',0,'FloatFormat','\t%.5g')
%
% license:
%     BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details
%
% -- this function is part of jsonlab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%

varname=inputname(2);
opt=varargin2struct(varargin{:});
rootisarray=0;
rootlevel=1;
forceroot=jsonopt('ForceRootName',0,opt);
if((isnumeric(obj) || islogical(obj) || ischar(obj) || isstruct(obj) || iscell(obj)) && isempty(rootname) && forceroot==0)
    rootisarray=1;
    rootlevel=0;
else
    if(isempty(rootname))
        rootname=varname;
    end
end
if((isstruct(obj) || iscell(obj))&& isempty(rootname) && forceroot)
    rootname='root';
end
json=obj2json(rootname,obj,rootlevel,opt);
if(rootisarray)
    json=sprintf('%s\n',json);
else
    json=sprintf('{\n%s\n}\n',json);
end

%%-------------------------------------------------------------------------
function txt=obj2json(name,item,level,varargin)

if(iscell(item))
    txt=cell2json(name,item,level,varargin{:});
elseif(isstruct(item))
    txt=struct2json(name,item,level,varargin{:});
elseif(ischar(item))
    txt=str2json(name,item,level,varargin{:});
else
    txt=mat2json(name,item,level,varargin{:});
end

%%-------------------------------------------------------------------------
function txt=cell2json(name,item,level,varargin)
txt='';
if(~iscell(item))
        error('input is not a cell');
end

dim=size(item);
len=numel(item); % let's handle 1D cell first
padding1=repmat(sprintf('\t'),1,level-1);
padding0=repmat(sprintf('\t'),1,level);
if(len>1) 
    if(~isempty(name))
        txt=sprintf('%s"%s": [\n',padding0, name); name=''; 
    else
        txt=sprintf('%s[\n',padding0); 
    end
end
for i=1:len
    txt=sprintf('%s%s%s',txt,padding1,obj2json(name,item{i},level+(len>1),varargin{:}));
    if(i<len) txt=sprintf('%s%s',txt,sprintf(',\n')); end
end
if(len>1) txt=sprintf('%s\n%s]',txt,padding0); end

%%-------------------------------------------------------------------------
function txt=struct2json(name,item,level,varargin)
txt='';
if(~isstruct(item))
	error('input is not a struct');
end
len=numel(item);
padding1=repmat(sprintf('\t'),1,level-1);
padding0=repmat(sprintf('\t'),1,level);
sep=',';

if(~isempty(name)) 
    if(len>1) txt=sprintf('%s"%s": [\n',padding0,name); end
else
    if(len>1) txt=sprintf('%s[\n',padding0); end
end
for e=1:len
  names = fieldnames(item(e));
  if(~isempty(name) && len==1)
        txt=sprintf('%s%s"%s": {\n',txt,repmat(sprintf('\t'),1,level+(len>1)), name); 
  else
        txt=sprintf('%s%s{\n',txt,repmat(sprintf('\t'),1,level+(len>1))); 
  end
  if(~isempty(names))
    for i=1:length(names)
	txt=sprintf('%s%s',txt,obj2json(names{i},getfield(item(e),...
             names{i}),level+1+(len>1),varargin{:}));
        if(i<length(names)) txt=sprintf('%s%s',txt,','); end
        txt=sprintf('%s%s',txt,sprintf('\n'));
    end
  end
  txt=sprintf('%s%s}',txt,repmat(sprintf('\t'),1,level+(len>1)));
  if(e==len) sep=''; end
  if(e<len) txt=sprintf('%s%s',txt,sprintf(',\n')); end
end
if(len>1) txt=sprintf('%s\n%s]',txt,padding0); end

%%-------------------------------------------------------------------------
function txt=str2json(name,item,level,varargin)
global isoct
txt='';
if(~ischar(item))
        error('input is not a string');
end
len=size(item,1);
sep=sprintf(',\n');

padding1=repmat(sprintf('\t'),1,level);
padding0=repmat(sprintf('\t'),1,level+1);

if(~isempty(name)) 
    if(len>1) txt=sprintf('%s"%s": [\n',padding1,name); end
else
    if(len>1) txt=sprintf('%s[\n',padding1); end
end
for e=1:len
    if(isoct)
        val=regexprep(item(e,:),'\\','\\');
        val=regexprep(val,'"','\"');
        val=regexprep(val,'^"','\"');
    else
        val=regexprep(item(e,:),'\\','\\\\');
        val=regexprep(val,'"','\\"');
        val=regexprep(val,'^"','\\"');
    end
    if(len==1)
        obj=['"' name '": ' '"',val,'"'];
	if(isempty(name)) obj=['"',val,'"']; end
        txt=sprintf('%s%s%s%s',txt,repmat(sprintf('\t'),1,level),obj);
    else
        txt=sprintf('%s%s%s%s',txt,repmat(sprintf('\t'),1,level+1),['"',val,'"']);
    end
    if(e==len) sep=''; end
    txt=sprintf('%s%s',txt,sep);
end
if(len>1) txt=sprintf('%s\n%s%s',txt,padding1,']'); end

%%-------------------------------------------------------------------------
function txt=mat2json(name,item,level,varargin)
if(~isnumeric(item) && ~islogical(item))
        error('input is not an array');
end

padding1=repmat(sprintf('\t'),1,level);
padding0=repmat(sprintf('\t'),1,level+1);

if(length(size(item))>2 || issparse(item) || ~isreal(item) || jsonopt('ArrayToStruct',0,varargin{:}))
    if(isempty(name))
    	txt=sprintf('%s{\n%s"_ArrayType_": "%s",\n%s"_ArraySize_": %s,\n',...
              padding1,padding0,class(item),padding0,regexprep(mat2str(size(item)),'\s+',',') );
    else
    	txt=sprintf('%s"%s": {\n%s"_ArrayType_": "%s",\n%s"_ArraySize_": %s,\n',...
              padding1,name,padding0,class(item),padding0,regexprep(mat2str(size(item)),'\s+',',') );
    end
else
    if(isempty(name))
    	txt=sprintf('%s%s',padding1,matdata2json(item,level+1,varargin{:}));
    else
        if(numel(item)==1 && jsonopt('NoRowBracket',1,varargin{:})==1)
            numtxt=regexprep(regexprep(matdata2json(item,level+1,varargin{:}),'^\[',''),']','');
           	txt=sprintf('%s"%s": %s',padding1,name,numtxt);
        else
    	    txt=sprintf('%s"%s": %s',padding1,name,matdata2json(item,level+1,varargin{:}));
        end
    end
    return;
end
dataformat='%s%s%s%s%s';

if(issparse(item))
    [ix,iy]=find(item);
    data=full(item(find(item)));
    if(~isreal(item))
       data=[real(data(:)),imag(data(:))];
       txt=sprintf(dataformat,txt,padding0,'"_ArrayIsComplex_": ','1', sprintf(',\n'));
    end
    txt=sprintf(dataformat,txt,padding0,'"_ArrayIsSparse_": ','1', sprintf(',\n'));
    if(find(size(item)==1))
        txt=sprintf(dataformat,txt,padding0,'"_ArrayData_": ',...
           matdata2json([ix,data],level+2,varargin{:}), sprintf('\n'));
    else
        txt=sprintf(dataformat,txt,padding0,'"_ArrayData_": ',...
           matdata2json([ix,iy,data],level+2,varargin{:}), sprintf('\n'));
    end
else
    if(isreal(item))
        txt=sprintf(dataformat,txt,padding0,'"_ArrayData_": ',...
            matdata2json(item(:)',level+2,varargin{:}), sprintf('\n'));
    else
        txt=sprintf(dataformat,txt,padding0,'"_ArrayIsComplex_": ','1', sprintf(',\n'));
        txt=sprintf(dataformat,txt,padding0,'"_ArrayData_": ',...
            matdata2json([real(item(:)) imag(item(:))],level+2,varargin{:}), sprintf('\n'));
    end
end
txt=sprintf('%s%s%s',txt,padding1,'}');

%%-------------------------------------------------------------------------
function txt=matdata2json(mat,level,varargin)
if(size(mat,1)==1)
    pre='';
    post='';
    level=level-1;
else
    pre=sprintf('[\n');
    post=sprintf('\n%s]',repmat(sprintf('\t'),1,level-1));
end
if(isempty(mat))
    txt='null';
    return;
end
floatformat=jsonopt('FloatFormat','%.10g',varargin{:});
formatstr=['[' repmat([floatformat ','],1,size(mat,2)-1) [floatformat sprintf('],\n')]];

if(nargin>=2 && size(mat,1)>1 && jsonopt('ArrayIndent',1,varargin{:})==1)
    formatstr=[repmat(sprintf('\t'),1,level) formatstr];
end
txt=sprintf(formatstr,mat');
txt(end-1:end)=[];
if(islogical(mat) && jsonopt('ParseLogical',0,varargin{:})==1)
   txt=regexprep(txt,'1','true');
   txt=regexprep(txt,'0','false');
end
%txt=regexprep(mat2str(mat),'\s+',',');
%txt=regexprep(txt,';',sprintf('],\n['));
% if(nargin>=2 && size(mat,1)>1)
%     txt=regexprep(txt,'\[',[repmat(sprintf('\t'),1,level) '[']);
% end
txt=[pre txt post];
if(any(isinf(mat(:))))
    txt=regexprep(txt,'([-+]*)Inf','"$1_Inf_"');
end
if(any(isnan(mat(:))))
    txt=regexprep(txt,'NaN','"_NaN_"');
end

%%-------------------------------------------------------------------------
function val=jsonopt(key,default,varargin)
val=default;
if(nargin<=2) return; end
opt=varargin{1};
if(isstruct(opt) && isfield(opt,key))
    val=getfield(opt,key);
end
