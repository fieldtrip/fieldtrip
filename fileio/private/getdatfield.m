function [datfield, dimord] = getdatfield(data)

% GETDATFIELD
%
% Use as
%   [datfield, dimord] = getdatfield(data)
% where the output arguments are cell-arrays.
%
% See also GETDIMORD, GETDIMSIZ

datfield  = fieldnames(data);

% these descriptive fields are cell-arrays and not treated as data
% the descriptive fields such as time and freq can be treated as data
xtrafield =  {'label' 'labelcmb'};
datfield  = setdiff(datfield, xtrafield);

xtrafield =  {'cfg' 'hdr' 'fsample' 'fsampleorig' 'grad' 'elec' 'opto' 'transform' 'dim' 'unit' 'coordsys' 'tri' 'tet' 'hex'};
datfield  = setdiff(datfield, xtrafield);

orgdim1   = datfield(~cellfun(@isempty, regexp(datfield, 'label$'))); % xxxlabel
datfield  = setdiff(datfield, orgdim1);
datfield  = datfield(:)';

orgdim1   = datfield(~cellfun(@isempty, regexp(datfield, 'dimord$'))); % xxxdimord
datfield  = setdiff(datfield, orgdim1);
datfield  = datfield(:)';

dimord = cell(size(datfield));
for i=1:length(datfield)
  dimord{i} = getdimord(data, datfield{i});
end
