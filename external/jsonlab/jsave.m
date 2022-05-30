function varargout=jsave(filename, varargin)
%
% jsave
%   or
% jsave(fname)
% varlist=jsave(fname,'param1',value1,'param2',value2,...)
%
% Store variables in a workspace to a JSON or binary JSON file
%
% authors:Qianqian Fang (q.fang <at> neu.edu)
% created on 2020/05/31
%
% input:
%      fname: (optional) output file name; if not given, save to 'jamdata.jamm'
%           if fname has a '.json' or '.jdt' suffix, a text-based
%           JSON/JData file will be created (slow); if the suffix is '.jamm' or
%           '.jdb', a Binary JData (https://github.com/NeuroJSON/bjdata/) file will be created.
%      opt: (optional) a struct to store parsing options, opt can be replaced by 
%           a list of ('param',value) pairs - the param string is equivallent
%           to a field in opt. opt can have the following 
%           fields (first in [.|.] is the default)
%
%           ws ['caller'|'base']: the name of the workspace in which the
%                         variables are to be saved
%           vars [{'var1','var2',...}]: cell array of variable names to be saved
%           matlab [0|1] if set to 1, use matlab's built-in jsonencode to
%                         store encoded data to a json file; output file
%                         must have a suffix of .jdt
%
%           all options for savebj/savejson (depends on file suffix)
%           can be used to adjust the output unless "'matlab',1" is used
%
% output:
%      varlist: a list of variables loaded
%
% examples:
%      jsave  % save all variables in the 'caller' workspace to jamdata.jamm
%      jsave('mydat.jamm','vars', {'v1','v2',...}) % save selected variables
%      jsave('mydat.jamm','compression','lzma')
%
% license:
%     BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details 
%
% -- this function is part of JSONLab toolbox (http://openjdata.org/jsonlab)
%

if(nargin==0)
    filename=[pwd filesep 'jamdata.jamm'];
end

opt=varargin2struct(varargin{:});

ws=jsonopt('ws','caller',opt);

allvar=evalin(ws,'whos');
varlist=jsonopt('vars',{allvar.name},opt);

[isfound, dontsave]=ismember(varlist,{allvar.name});
if(any(isfound==0))
    error('specified variable is not found');
end

header=struct;
body=struct;

metadata=struct('CreateDate',datestr(now,29),...
                'CreateTime',datestr(now,'hh:mm:ss'),...
                'OriginalName',filename,...
                'BJDataVersion','v1_draft-2');

vers=ver('MATLAB');
if(isempty(vers))
    vers=ver('Octave');
    [verstr, releasedate]=version;
    vers.Release=verstr;
    vers.Date=releasedate;
end

metadata.CreatorApp=vers.Name;
metadata.CreatorVersion=vers.Version;
metadata.CreatorRelease=vers.Release;
metadata.ReleaseDate=vers.Date;
metadata.FormatVersion=1;
metadata.Parameters=opt;

header.(encodevarname('_DataInfo_'))=metadata;

for i=1:length(varlist)
    header.(varlist{i})=allvar(dontsave(i));
    body.(varlist{i})=evalin(ws,varlist{i});
end

savefun=@savebj;
loadfun=@loadbj;
if(regexp(filename,'\.[jJ][sS][oO][nN]$'))
    savefun=@savejson;
    loadfun=@loadjson;
elseif(regexp(filename,'\.[jJ][dD][tT]$'))
    savefun=@savejson;
    loadfun=@loadjson;
elseif(regexp(filename,'\.[mM][sS][gG][pP][kK]$'))
    savefun=@savemsgpack;
end

if(nargout==1)
    varargout{1}=header;
end

if(jsonopt('matlab',0,opt) && exist('jsonencode','builtin'))
    if(isempty(regexp(filename,'\.[jJ][sS][oO][nN]$', 'once')))
        filename=regexprep(filename,'\.[a-zA-Z]+$','.jdt');
    end
    output.WorkspaceHeader=jdataencode(header,'prefix','x','base64',1,varargin{:});
    headerjson=jsonencode(output);
    clear output;

    output.WorkspaceData=jdataencode(body,'AnnotateArray',1,'base64',1,...
       'Compression','zlib','UseArrayZipSize',1,'MapAsStruct',1,...
       'prefix','x',varargin{:});
    bodyjson=jsonencode(output);
    clear output;

    fid=fopen(filename,jsonopt('writemode','w',opt));
    fwrite(fid,headerjson);
    fwrite(fid,sprintf('\n\n\n'));
    fwrite(fid,bodyjson);
    fclose(fid);
elseif(jsonopt('usemmap',0,opt)==1)
    bodyjson=savefun('WorkspaceData',body,...
        'compression','zlib','keeptype',1,'array2struct',1,varargin{:});
    header.(encodevarname('_MMap_'))=loadfun(bodyjson,'mmaponly',1,varargin{:});
    savefun('WorkspaceHeader',header,'filename',filename,varargin{:});
    fid=fopen(filename,'ab+');
    fwrite(fid,bodyjson);
    fclose(fid);
else
    savefun('WorkspaceHeader',header,'filename',filename,varargin{:});
    savefun('WorkspaceData',body,'filename',filename,'append',1,...
        'compression','zlib','keeptype',1,'array2struct',1,varargin{:});
end
