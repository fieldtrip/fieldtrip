function jsn = jsnirfcreate(varargin)
%
%    jsn=jsnirfcreate
%       or
%    jsn=jsnirfcreate(option)
%    jsn=jsnirfcreate('Format',format,'Param1',value1, 'Param2',value2,...)
%
%    Create an empty JSNIRF data structure defined in the JSNIRF
%    specification: https://github.com/NeuroJSON/jsnirf or a SNIRF data structure
%    based on https://github.com/fNIRS/snirf
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        option (optional): option can be ignored. If it is a string with a
%             value 'snirf', this creates a default SNIRF data structure;
%             otherwise, a JSNIRF data structure is created.
%        format: same as option.
%        param/value:   a list of name/value pairs specify
%             additional subfields to be stored under the /nirs object.
%
%    output:
%        jsn: a default SNIRF or JSNIRF data structure.
%
%    example:
%        jsn=jsnirfcreate('data',mydata,'aux',myauxdata,'comment','test');
%
%    this file is part of JSNIRF specification: https://github.com/NeuroJSON/jsnirf
%
%    License: GPLv3 or Apache 2.0, see https://github.com/NeuroJSON/jsnirf for details
%

% define empty SNIRF data structure with all required fields

defaultmeta = struct('SubjectID', 'default', 'MeasurementDate', datestr(now, 29), ...
                     'MeasurementTime', datestr(now, 'hh:mm:ss'), 'LengthUnit', 'mm', ...
                     'TimeUnit', 's', 'FrequencyUnit', 'Hz');
defaultsrcmap = struct('sourceIndex', [], 'detectorIndex', [], ...
                       'wavelengthIndex', [], 'dataType', 1, 'dataTypeIndex', 1);
defaultdata = struct('dataTimeSeries', [], 'time', [], 'measurementList', defaultsrcmap);
defaultaux = struct('name', '', 'dataTimeSeries', [], 'time', [], 'timeOffset', 0);
defaultstim = struct('name', '', 'data', []);
defaultprobe = struct('wavelengths', [], 'sourcePos2D', [], 'detectorPos2D', []);

nirsdata = struct('metaDataTags', defaultmeta, ...
                  'data', defaultdata, ...
                  'aux', defaultaux, ...
                  'stim', defaultstim, ...
                  'probe', defaultprobe);

% read user specified data fields - will validate format in future updates

if (nargin > 1 && bitand(nargin, 1) == 0)
    for i = 1:nargin * 0.5
        key = varargin{2 * i - 1};
        if (strcmpi(key, 'format'))
            key = 'format';
        end
        nirsdata.(key) = varargin{2 * i};
    end
end

jsn = struct();

% return either a SNIRF data structure, or JSNIRF data (enclosed in SNIRFData tag)

if ((nargin == 1 && strcmpi(varargin{1}, 'snirf')) || ...
    (isfield(nirsdata, 'format') && strcmpi(nirsdata.format, 'snirf')))
    if (isfield(nirsdata, 'format'))
        nirsdata = rmfield(nirsdata, 'format');
    end
    jsn = struct('formatVersion', '1.0', 'nirs', nirsdata);
else
    nirsdata.formatVersion = '1.0';
    len = length(fieldnames(nirsdata));
    nirsdata = orderfields(nirsdata, [len, 1:len - 1]);
    jsn = struct('SNIRFData', nirsdata);
end
