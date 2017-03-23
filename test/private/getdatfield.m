function [datfield, dimord] = getdatfield(data)

% GETDATFIELD
%
% Use as
%   [datfield, dimord] = getdatfield(data)
% where the output arguments are cell-arrays.
%
% See also GETDIMORD, GETDIMSIZ

datfield  = fieldnames(data);
datfield  = datfield(:)'; % it should be a row-vector throughout the subsequent code

% these descriptive fields are cell-arrays and not treated as data
% the descriptive fields such as time and freq can be treated as data
xtrafield =  {'label' 'labelcmb'};
datfield  = setdiff(datfield, xtrafield);

xtrafield =  {'cfg' 'hdr' 'fsample' 'fsampleorig' 'grad' 'elec' 'opto' 'transform' 'transformorig' 'dim' 'unit' 'coordsys' 'tri' 'tet' 'hex'};
datfield  = setdiff(datfield, xtrafield);

% find substructure fields
structfield = {};
for i=1:numel(datfield)
  if isstruct(data.(datfield{i}))
    structfield(end+1) = datfield(i);
  end
end
datfield  = setdiff(datfield, structfield);

% find homogenous transformation matrices, i.e. a 4x4 matrix that is named xxx2yyy
transformfield = {};
for i=1:numel(datfield)
  if ~isempty(regexp(datfield, '^[a-zA-Z]*2[a-zA-Z]*$', 'once')) && isequal(size(data.(datfield{i})), [4 4])
    transformfield(end+1) = datfield(i);
  end
end
datfield  = setdiff(datfield, transformfield);

xxxlabel  = datfield(~cellfun(@isempty, regexp(datfield, 'label$'))); % xxxlabel
datfield  = setdiff(datfield, xxxlabel);

xxxdimord = datfield(~cellfun(@isempty, regexp(datfield, 'dimord$'))); % xxxdimord
datfield  = setdiff(datfield, xxxdimord);

dimord = cell(size(datfield));
for i=1:length(datfield)
  dimord{i} = getdimord(data, datfield{i});
end

