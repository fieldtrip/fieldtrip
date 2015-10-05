function dimsiz = getdimsiz(data, field)

% GETDIMSIZ
%
% Use as
%   dimsiz = getdimsiz(data, field)
%
% If the length of the vector that is returned is smaller than the
% number of dimensions that you would expect from GETDIMORD, you
% should assume that it has trailing singleton dimensions.
%
% Example use
%   dimord = getdimord(datastructure, fieldname);
%   dimtok = tokenize(dimord, '_');
%   dimsiz = getdimsiz(datastructure, fieldname);
%   dimsiz(end+1:length(dimtok)) = 1; % there can be additional trailing singleton dimensions
%
% See also GETDIMORD, GETDATFIELD

if ~isfield(data, field) && isfield(data, 'avg') && isfield(data.avg, field)
  field = ['avg.' field];
elseif ~isfield(data, field) && isfield(data, 'trial') && isfield(data.trial, field)
  field = ['trial.' field];
elseif ~isfield(data, field)
  error('field "%s" not present in data', field);
end

if strncmp(field, 'avg.', 4)
  prefix = [];
  field = field(5:end); % strip the avg
  data.(field) = data.avg.(field); % move the avg into the main structure
  data = rmfield(data, 'avg');
elseif strncmp(field, 'trial.', 6)
  prefix = numel(data.trial);
  field = field(7:end); % strip the trial
  data.(field) = data.trial(1).(field); % move the first trial into the main structure
  data = rmfield(data, 'trial');
else
  prefix = [];
end

dimsiz = cellmatsize(data.(field));

% add nrpt in case of source.trial
dimsiz = [prefix dimsiz];

end % main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to determine the size of data representations like {pos}_ori_time
% FIXME this will fail for {xxx_yyy}_zzz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function siz = cellmatsize(x)
if iscell(x)
  if isempty(x)
    siz = 0;
    return % nothing else to do
  elseif isvector(x)
    cellsize = numel(x);          % the number of elements in the cell-array
  else
    cellsize = size(x);
    x = x(:); % convert to vector for further size detection
  end
  [dum, indx] = max(cellfun(@numel, x));
  matsize = size(x{indx});      % the size of the content of the cell-array
  siz = [cellsize matsize];     % concatenate the two
else
  siz = size(x);
end
end % function cellmatsize
