function [boolval, list] = checkchan(varargin)

% last input is always the required string
required = varargin{end};
varargin = varargin(1:end-1);

Ndata = numel(varargin);

if strcmp(required, 'unique')
  Nchan = zeros(1,Ndata);
  list  = cell(0,1);
  for i=1:Ndata
    Nchan(i) = numel(varargin{i}.label);
    list     = [list;varargin{i}.label(:)];
  end
  boolval = numel(unique(list))==numel(list);
  
elseif strcmp(required, 'identical')
  
  % determine whether channels are equal across the inputs
  % channel order is not yet handled here
  for k = 2:Ndata
    matched = numel(match_str(varargin{1}.label, varargin{k}.label));
    if matched ~= numel(varargin{1}.label) || matched ~= numel(varargin{k}.label)
      boolval = 0;
      return;
    end
  end
  
  boolval = 1;
end

