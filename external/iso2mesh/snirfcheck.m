function data = snirfcheck(root, varargin)
%
%    data=snirfcheck(root)
%       or
%    data=snirfcheck(root,type)
%    data=snirfcheck(root,{'nameA','nameB',...})
%
%    Processing an HDF5 based SNIRF data and group indexed datasets into a
%    cell array
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        root: the raw input snirf data structure (loaded from loadh5.m)
%        type: if type is set as a cell array of strings, it restrict the
%              grouping only to the subset of field names in this list;
%              if type is a string as 'snirf', it is the same as setting
%              type as {'aux','data','nirs','stim','measurementList'}.
%
%    output:
%        data: a reorganized matlab structure. Each SNIRF data chunk is
%              enclosed inside a 'SNIRFData' subfield or cell array.
%
%    example:
%        rawdata=loadh5('mydata.snirf');
%        data=snirfdecode(rawdata);
%
%    this file is part of JSNIRF specification: https://github.com/NeuroJSON/jsnirf
%
%    License: Apache 2.0, see https://github.com/NeuroJSON/jsnirf for details
%

if (nargin < 1)
    help snirfdecode;
    return
end

data = regrouph5(root, varargin{:});

issnirf = 1;

if (~isempty(varargin))
    if (ischar(varargin{1}) && strcmpi(varargin{1}, 'jsnirf'))
        issnirf = 0;
    end
end

if (issnirf == 0 && isfield(data, 'nirs') && isfield(data, 'formatVersion') && ~isfield(data, 'SNIRFData'))
    data.SNIRFData = data.nirs;
    if (isfield(data.SNIRFData, 'data') && isfield(data.SNIRFData.data, 'measurementList'))
        data.SNIRFData.data.measurementList = aos2soa(data.nirs.data.measurementList);
    end
    if (iscell(data.nirs))
        for i = 1:length(data.nirs)
            data.SNIRFData{i}.formatVersion = data.formatVersion;
            len = length(fieldnames(data.SNIRFData{i}));
            data.SNIRFData{i} = orderfields(data.SNIRFData{i}, [len, 1:len - 1]);
        end
    else
        data.SNIRFData.formatVersion = data.formatVersion;
        len = length(fieldnames(data.SNIRFData));
        data.SNIRFData = orderfields(data.SNIRFData, [len, 1:len - 1]);
    end
    data = rmfield(data, {'nirs', 'formatVersion'});
end
