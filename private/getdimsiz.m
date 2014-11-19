function dimsiz = getdimsiz(data, field)

% GETDIMSIZ
%
% Use as
%   dimsiz = getdimsiz(data, field)
%
% See also GETDIMORD

if strncmp(field, 'avg.', 4)
  field = field(5:end); % strip the avg
  data.(field) = data.avg.(field);
end

if ~isfield(data, field)
  error('field "%s" not present in data', field);
end

dimsiz = cellmatsize(data.(field));

end % main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to determine the size of data representations like {pos}_ori_time
% FIXME this will fail for {xxx_yyy}_zzz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function siz = cellmatsize(x)
if iscell(x)
  cellsize = numel(x);          % the number of elements in the cell-array
  [dum, indx] = max(cellfun(@numel, x));
  matsize = size(x{indx});      % the size of the content of the cell-array
  siz = [cellsize matsize];     % concatenate the two
else
  siz = size(x);
end
end % function cellmatsize