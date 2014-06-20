function tmpcfg = maketmpcfg(cfg, fields)

% MAKETMPCFG copies all specified fields that are present in the input cfg
% to an empty structure, which is returned as tmpcfg.

if ~iscell(fields) && ischar(fields)
  fields = {fields};
elseif ~iscell(fields)
  error('fields input argument must be a cell array of strings or a single string');
end

tmpcfg = [];

for k = 1:numel(fields)
  if isfield(cfg, fields{k})
    tmpcfg.(fields{k}) = cfg.(fields{k});
  end
end

end