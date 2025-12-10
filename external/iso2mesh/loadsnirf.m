function data = loadsnirf(fname, varargin)
%
%    data=loadsnirf(fname)
%       or
%    jnirs=loadsnirf(fname, 'Param1',value1, 'Param2',value2,...)
%
%    Load an HDF5 based SNIRF file, and optionally convert it to a JSON
%    file based on the JSNIRF specification:
%    https://github.com/NeuroJSON/jsnirf
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        fname: the input snirf data file name (HDF5 based)
%
%    output:
%        data: a MATLAB structure with the grouped data fields
%
%    dependency:
%        - the loadh5/regrouph5 functions are provided by the eazyh5
%          toolbox at http://github.com/fangq/eazyh5
%        - the varargin2struct and jsonopt functions are provided by the JSONLab
%          toolbox at http://github.com/NeuroJSON/jsonlab
%        - if data compression is specified by 'compression','zlib' param/value
%          pairs, ZMat toolbox will be needed, http://github.com/fangq/zmat
%
%    example:
%        data=loadsnirf('test.snirf');
%
%    this file is part of JSNIRF specification: https://github.com/NeuroJSON/jsnirf
%
%    License: GPLv3 or Apache 2.0, see https://github.com/NeuroJSON/jsnirf for details
%

if (nargin == 0 || ~ischar(fname))
    error('you must provide a file name');
end

data = loadh5(fname, 'stringarray', 1);

opt = struct;

if (length(varargin) == 1)
    data = snirfdecode(data, varargin{:});
elseif (length(varargin) >= 2)
    opt = varargin2struct(varargin{:});
    data = snirfdecode(data, varargin{:});
else
    data = snirfdecode(data);
end

outfile = jsonopt('FileName', '', opt);
if (~isempty(outfile))
    if (regexp(outfile, '\.[Bb][Nn][Ii][Rr][Ss]$'))
        savebj('SNIRFData', data, 'FileName', outfile, opt);
    elseif (~isempty(regexp(outfile, '\.[Jj][Nn][Ii][Rr][Ss]$', 'once')) || ~isempty(regexp(outfile, '\.[Jj][Ss][Oo][Nn]$', 'once')))
        savejson('SNIRFData', data, 'FileName', outfile, opt);
    elseif (regexp(outfile, '\.[Mm][Aa][Tt]$'))
        save(outfile, 'data');
    else
        error('only support .jnirs,.bnirs and .mat files');
    end
end
