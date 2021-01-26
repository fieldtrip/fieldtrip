function val = HDF5_PostProcessing(val, val0, options)

if ~exist('val0','var')
    val0 = [];
end
if ~exist('options','var')
    options = '';
end
val = HDF5_Transpose(val, options);

% Convert muti-row char array to cell string array
if ischar(val) && ~iscell(val0)
    val = convertHDF5StrToMatlabStr(val);
elseif ischar(val)
    val = convertHDF5StrToMatlabStr(val, 'cell');
end
