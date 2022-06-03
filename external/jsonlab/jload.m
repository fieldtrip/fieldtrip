function varargout=jload(filename, varargin)
%
% jload
%   or
% jload(fname)
% varlist=jload(fname)
% [varlist, header]=jload(fname)
% varlist=jload(fname,'param1',value1,'param2',value2,...)
%
% Load variables from a JSON or binary JSON file to a workspace
%
% authors:Qianqian Fang (q.fang <at> neu.edu)
% created on 2020/05/31
%
% input:
%      fname: (optional) input file name; if not given, load 'jamdata.jamm'
%           if fname has a '.json' or '.jdt' suffix, a text-based
%           JSON/JData file will be expected; if the suffix is '.jamm' or
%           '.jdb', a Binary JData file will be expected.
%      opt: (optional) a struct to store parsing options, opt can be replaced by 
%           a list of ('param',value) pairs - the param string is equivallent
%           to a field in opt. opt can have the following 
%           fields (first in [.|.] is the default)
%
%           ws ['caller'|'base']: the name of the workspace in which the
%                         variables are to be saved
%           vars [{'var1','var2',...}]: list of variables to be saved
%           header [0|1]: if set to 1, return the metadata of the variables 
%                         stored in the file
%           matlab [0|1] if set to 1, use matlab's built-in jsondecode to
%                         parse the json file and then decode the output by
%                         jdatadecode; input file must have a suffix of .jdt
%
%           all options for loadbj/loadjson (depends on file suffix)
%           can be used to adjust the parsing options
%
% output:
%      varlist: a struct with each subfield a variable stored in the file,
%               if output is ignored, the variables will be loaded to the
%               workspace specified by the 'ws' option, which by default
%               load the variables to the current workspace ('caller')
%
% examples:
%      jload  % load all variables in jamdata.jamm to the 'caller' workspace 
%      jload mydat.jamm
%      jload('mydat.jamm','vars', {'v1','v2',...}) % load selected variables
%      varlist=jload('mydat.jamm','simplifycell',1)
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

loadfun=@loadbj;
if(regexp(filename,'\.[jJ][sS][oO][nN]$'))
    loadfun=@loadjson;
elseif(regexp(filename,'\.[jJ][dD][tT]$'))
    loadfun=@loadjson;
elseif(regexp(filename,'\.[mM][sS][gG][pP][kK]$'))
    loadfun=@loadmsgpack;
end

if(jsonopt('matlab',0,opt) && exist('jsonencode','builtin'))
    jsonstr=fileread(filename);
    pos=regexp(jsonstr,'}\n\n\n{"WorkspaceData":','once');
    if(isempty(pos))
        error('the json file is not generated using matlab''s jsonencode');
    end
    header=jsondecode(jsonstr(1:pos+1));
else
    try
        header=loadfun(filename,'ObjectID',1, 'MaxBuffer', 65536, varargin{:});
    catch
        header=loadfun(filename,'ObjectID',1, varargin{:});
    end
end

allvar=fieldnames(header.WorkspaceHeader);

varlist=jsonopt('vars',allvar,opt);
varlist(ismember(varlist,encodevarname('_DataInfo_')))=[];

isfound=ismember(varlist,allvar);
if(any(isfound==0))
    error('specified variable is not found');
end

if(jsonopt('matlab',0,opt) && exist('jsonencode','builtin'))
    body=jdatadecode(jsondecode(jsonstr(pos+4:end)));
else
    body=loadfun(filename,'ObjectID',2, varargin{:});
end

if(nargout==0)
    for i=1:length(varlist)
        assignin(ws, varlist{i}, body.WorkspaceData.(varlist{i}));
    end
else
    varargout{1}=rmfield(body.WorkspaceData,setdiff(fieldnames(body.WorkspaceData),varlist));
    if(nargout>1)
        varargout{2}=header;
    end
end
