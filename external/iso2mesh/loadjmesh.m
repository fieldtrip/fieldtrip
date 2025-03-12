function jmsh = loadjmesh(filename, varargin)
%
%    jmsh=loadjmesh(inputfile)
%       or
%    jmsh=loadjmesh(inputfile, 'Param1',value1, 'Param2',value2,...)
%
%    Load a text or binary JMesh file with format defined in JMesh
%    specification: https://github.com/NeuroJSON/jmesh
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        inputfile: the output file name to the JMesh
%                *.bmsh for binary JMesh file
%                *.jmsh for text JMesh file
%        options: (optional) if loading from a .bmsh file, please see the options for
%               loadbj.m (part of JSONLab); if loading from a .jmsh, please see the
%               supported options for loadjson.m (part of JSONLab).
%
%    output:
%        jmsh: a structure (array) or cell (array). The data structure can
%            be completely generic or auxilary data without any JMesh
%            constructs. However, if a JMesh object is included, it shall
%            contain the constructs defined in the JMesh specification
%
%    example:
%        [no,fc,el]=meshasphere([0 0 0],50,3,10);
%        savejmesh(no,fc,el,'box.jmsh');
%        newmesh=loadjmesh('box.jmsh');
%
%    this file is part of JMesh specification: https://github.com/NeuroJSON/jmesh
%
%    License: Apache 2.0, see https://github.com/NeuroJSON/jmesh for details
%

if (nargin < 1)
    error('you must provide data and output file name');
end

if (~exist('loadbj', 'file'))
    error('you must first install JSONLab from http://github.com/NeuroJSON/jsonlab/');
end

if (regexp(filename, '\.jmsh$'))
    jmsh = loadjson(filename, varargin{:});
elseif (regexp(filename, '\.bmsh$'))
    jmsh = loadbj(filename, varargin{:});
else
    error('file suffix must be .jmsh for text JMesh, .bmsh for binary JMesh');
end
