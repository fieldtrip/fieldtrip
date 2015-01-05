function dimsiz = getdimsiz(data, field)

% GETDIMSIZ
%
% Use as
%   dimsiz = getdimsiz(data, field)
%
% See also GETDIMORD

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
  data.(field) = data.avg.(field); % copy the avg into the main structure
  data = rmfield(data, 'avg');
elseif strncmp(field, 'trial.', 6)
  prefix = numel(data.trial);
  field = field(7:end); % strip the trial
  data.(field) = data.trial(1).(field); % copy the first trial into the main structure
  data = rmfield(data, 'trial');
else
  prefix = [];
end

dimsiz = cellmatsize(data.(field));

% figure out whether there are additional trailing singleton dimensions
try
  % don't call getdimord in case getdimord is calling us, we don't want to recurse
  s = dbstack;
  if ~strcmp(s(2).name, 'getdimord')
    dimord = getdimord(data, field);
    dimtok = tokenize(dimord, '_');
    for i=(length(dimsiz)+1):length(dimtok)
      dimsiz(i) = 1;
    end
  end
end % try

% add nrpt in case of source.trial
dimsiz = [prefix dimsiz];

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