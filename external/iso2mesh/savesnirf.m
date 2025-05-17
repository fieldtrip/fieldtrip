function savesnirf(data, outfile, varargin)
%
%    savesnirf(snirfdata, fname)
%       or
%    savesnirf(snirfdata, fname, 'Param1',value1, 'Param2',value2,...)
%
%    Load an HDF5 based SNIRF file, and optionally convert it to a JSON
%    file based on the JSNIRF specification:
%    https://github.com/NeuroJSON/jsnirf
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        snirfdata: a raw SNIRF data, preprocessed SNIRF data or JSNIRF
%             data (root object must be SNIRFData)
%        fname: the output SNIRF (.snirf) or JSNIRF data file name (.jnirs, .bnirs)
%
%    output:
%        data: a MATLAB structure with the grouped data fields
%
%    example:
%        data=loadsnirf('test.snirf');
%        savesnirf(data,'newfile.snirf');
%
%    this file is part of JSNIRF specification: https://github.com/NeuroJSON/jsnirf
%
%    License: GPLv3 or Apache 2.0, see https://github.com/NeuroJSON/jsnirf for details
%

if (nargin < 2 || ~ischar(outfile))
    error('you must provide data and a file name');
end

opt = varargin2struct(varargin{:});
if (~isfield(opt, 'root'))
    opt.rootname = '';
end
if (~isfield(opt, 'variablelengthstring'))
    opt.variablelengthstring = 1;
end
if (~isfield(opt, 'rowas1d'))
    opt.rowas1d = 1;
end

if (isfield(data, 'SNIRFData'))
    data.nirs = data.SNIRFData;
    data.formatVersion = data.SNIRFData.formatVersion;
    data.nirs = rmfield(data.nirs, 'formatVersion');
    data = rmfield(data, 'SNIRFData');
end

if (~isempty(outfile))
    if (~isempty(regexp(outfile, '\.[Hh]5$', 'once')))
        saveh5(data, outfile, opt);
    elseif (~isempty(regexp(outfile, '\.[Ss][Nn][Ii][Rr][Ff]$', 'once')))
        if (isfield(data.nirs.data, 'measurementList'))
            forceint = {'sourceIndex', 'detectorIndex', 'wavelengthIndex', ...
                        'dataType', 'dataTypeIndex', 'moduleIndex', ...
                        'sourceModuleIndex', 'detectorModuleIndex'};
            for i = 1:length(forceint)
                if (isfield(data.nirs.data.measurementList, forceint{i}))
                    if (iscell(data.nirs.data.measurementList.(forceint{i})))
                        data.nirs.data.measurementList.(forceint{i}) = cell2mat(data.nirs.data.measurementList.(forceint{i}));
                    end
                    data.nirs.data.measurementList.(forceint{i}) = int32(data.nirs.data.measurementList.(forceint{i}));
                end
            end
            if (length(data.nirs.data.measurementList) == 1 && ...
                length(data.nirs.data.measurementList.sourceIndex) > 1)
                data.nirs.data.measurementList = soa2aos(data.nirs.data.measurementList);
            end
        end
        if (opt.rowas1d)
            force1d.probe = {'wavelengths', 'wavelengthsEmission', 'frequencies', ...
                             'timeDelays', 'timeDelayWidths', 'momentOrders', 'correlationTimeDelays', ...
                             'correlationTimeDelayWidths'};
            force1d.data = {'time'};
            force1d.aux = {'time', 'timeOffset'};
            fields = fieldnames(force1d);
            for i = 1:length(fields)
                for j = 1:length(force1d.(fields{i}))
                    if (isfield(data.nirs.(fields{i}), force1d.(fields{i}){j}))
                        if (iscell(data.nirs.(fields{i}).(force1d.(fields{i}){j})))
                            data.nirs.(fields{i}).(force1d.(fields{i}){j}) = cell2mat(data.nirs.(fields{i}).(force1d.(fields{i}){j}));
                        end
                        data.nirs.(fields{i}).(force1d.(fields{i}){j}) = timeseries(data.nirs.(fields{i}).(force1d.(fields{i}){j})(:).');
                    end
                end
            end
        end
        if (isfield(data.nirs, 'probe'))
            forcestrarray.probe = {'sourceLabels', 'detectorLabels', 'landmarkLabels'};
            forcestrarray.stim = {'dataLabels'};
            fields = fieldnames(forcestrarray);
            for i = 1:length(fields)
                for j = 1:length(forcestrarray.(fields{i}))
                    if (isfield(data.nirs.(fields{i}), forcestrarray.(fields{i}){j}))
                        if (iscell(data.nirs.(fields{i}).(forcestrarray.(fields{i}){j})))
                            data.nirs.(fields{i}).(forcestrarray.(fields{i}){j}) = cell2mat(data.nirs.(fields{i}).(forcestrarray.(fields{i}){j}));
                        end
                        data.nirs.(fields{i}).(forcestrarray.(fields{i}){j}) = timeseries(string(data.nirs.(fields{i}).(forcestrarray.(fields{i}){j})(:).'));
                    end
                end
            end
        end
        if (~isempty(regexp(data.formatVersion, '^1\.', 'once')))
            if (length(data.nirs.data.measurementList) == 1 && length(data.nirs.data.measurementList.sourceIndex) > 1)
                data.nirs.data.measurementList = soa2aos(data.nirs.data.measurementList);
            end
        end
        data.nirs.data = forceindex(data.nirs.data, 'measurementList');
        data.nirs = forceindex(data.nirs, 'data');
        data.nirs = forceindex(data.nirs, 'stim');
        data.nirs = forceindex(data.nirs, 'aux');
        saveh5(data, outfile, opt);
    elseif (~isempty(regexp(outfile, '\.[Jj][Nn][Ii][Rr][Ss]$', 'once')) || ~isempty(regexp(outfile, '\.[Jj][Ss][Oo][Nn]$', 'once')))
        savejson('SNIRFData', data, 'FileName', outfile, opt);
    elseif (regexp(outfile, '\.[Mm][Aa][Tt]$'))
        save(outfile, 'data');
    elseif (regexp(outfile, '\.[Bb][Nn][Ii][Rr][Ss]$'))
        savebj('SNIRFData', data, 'FileName', outfile, opt);
    else
        error('only support .snirf, .h5, .jnirs, .bnirs and .mat files');
    end
end

% force adding index 1 to the group name for singular struct and cell
function newroot = forceindex(root, name)
newroot = root;
fields = fieldnames(newroot);
idx = find(ismember(fields, name));
if (~isempty(idx) && length(newroot.(name)) == 1)
    newroot.(sprintf('%s1', name)) = newroot.(name);
    newroot = rmfield(newroot, name);
    fields{idx(1)} = sprintf('%s1', name);
    newroot = orderfields(newroot, fields);
end
