function jnirs = loadjsnirf(filename, varargin)
%
%    jnirs=loadjsnirf(inputfile)
%       or
%    jnirs=loadjsnirf(inputfile, 'Param1',value1, 'Param2',value2,...)
%
%    Load a text (.jnirs or .json) or binary (.bnirs) based JSNIRF
%    file defined in the JSNIRF specification:
%    https://github.com/NeuroJSON/jsnirf or a .snirf/.h5 SNIRF data defined in
%    the SNIRF specification https://github.com/fNIRS/snirf
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        inputfile: the output file name to the JSNIRF or SNIRF file
%                *.bnirs for binary JSNIRF file
%                *.jnirs for text JSNIRF file
%                *.snirf for HDF5/SNITRF files
%        options: (optional) if loading from a .bnii file, please see the options for
%               loadbj.m (part of JSONLab); if loading from a .jnirs, please see the
%               supported options for loadjson.m (part of JSONLab).
%
%    output:
%        jnirs: a structure (array) or cell (array). The data structure can
%            be completely generic or auxilary data without any JSNIRF
%            constructs. However, if a JSNIRF object is included, it shall
%            contain the below subfields (can appear within any depth of the
%            structure)
%
%    example:
%        newjnirs=loadjsnirf('subject1.jnirs');
%
%    this file is part of JSNIRF specification: https://github.com/NeuroJSON/jsnirf
%
%    License: GPLv3 or Apache 2.0, see https://github.com/NeuroJSON/jsnirf for details
%

if (nargin < 1)
    error('you must provide data and output file name');
end

if (~exist('savejson', 'file'))
    error('you must first install JSONLab from http://github.com/NeuroJSON/jsonlab/');
end

if (~isempty(regexp(filename, '\.[Ss][Nn][Ii][Rr][Ff]$', 'once')) || ~isempty(regexp(filename, '\.[Hh]5$', 'once')))
    jnirs = loadsnirf(filename, varargin);
elseif (regexp(filename, '\.[Jj][Nn][Ii][Rr][Ss]$'))
    jnirs = loadjson(filename, varargin{:});
elseif (regexp(filename, '\.[Bb][Nn][Ii][Rr][Ss]$'))
    jnirs = loadbj(filename, varargin{:});
else
    error('file suffix must be .snirf for SNIRF/HDF5, .jnirs for text JSNIRF, .bnirs for binary JSNIRF files');
end
