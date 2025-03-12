function data = savejsnirf(jnirs, filename, varargin)
%
%    savejsnirf(jnirs, outputfile)
%       or
%    savejsnirf(jnirs, outputfile, 'Param1',value1, 'Param2',value2,...)
%
%    Save an in-memory JSNIRF structure into a JSNIRF file with format
%    defined in JSNIRF specification: https://github.com/NeuroJSON/jsnirf
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        jnirs: a structure (array) or cell (array). The data structure can
%            be completely generic or auxilary data without any JSNIRF
%            constructs. However, if a JSNIRF object is included, it shall
%            contain the below subfields (can appear within any depth of the
%            structure)
%                jnirs.SNIRFData - the main image data array
%        outputfile: the output file name to the JSNIRF file
%                *.bnirs for binary JSNIRF file
%                *.jnirs for text JSNIRF file
%                *.snirf or *.h5 for HDF5-based SNIRF file
%        options: (optional) if saving to a .bnirs file, please see the options for
%               savebj.m (part of JSONLab); if saving to .jnirs, please see the
%               supported options for savejson.m (part of JSONLab).
%
%    dependency:
%        - the savejson and savebj functions are provided by the JSONLab
%          toolbox at http://github.com/NeuroJSON/jsonlab
%        - if data compression is specified by 'compression','zlib' param/value
%          pairs, ZMat toolbox will be needed, http://github.com/fangq/zmat
%        - the saveh5 function is provided by the eazyh5 toolbox at
%          http://github.com/fangq/eazyh5
%
%    example:
%        jnirs=jsnirfcreate('aux',struct('name','pO2','dataTimeSeries',1:10,'time',1:10));
%        savejsnirf(jnirs, 'test.jnirs');
%        savejsnirf(jnirs, 'test.bnirs','compression','zlib');
%
%    this file is part of JSNIRF specification: https://github.com/NeuroJSON/jsnirf
%
%    License: GPLv3 or Apache 2.0, see https://github.com/NeuroJSON/jsnirf for details
%

if (nargin < 2)
    error('you must provide data and output file name');
end

if (~exist('savejson', 'file'))
    error('you must first install JSONLab from http://github.com/NeuroJSON/jsonlab/');
end

data = '';

if (regexp(filename, '\.[Jj][Nn][Ii][Rr][Ss]$'))
    data = savejson('', jnirs, 'filename', filename, varargin{:});
elseif (regexp(filename, '\.[Bb][Nn][Ii][Rr][Ss]$'))
    data = savebj('', jnirs, 'filename', filename, varargin{:});
elseif (~isempty(regexp(filename, '\.[Ss][Nn][Ii][Rr][Ff]$', 'once')) || ~isempty(regexp(outfile, '\.[Hh]5$', 'once')))
    saveh5(jnirs, filename, varargin{:});
else
    error('file suffix must be .jnirs for text JSNIRF or .bnirs for binary JSNIRF');
end
