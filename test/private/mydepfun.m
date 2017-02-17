function [outlist, depmat] = mydepfun(inlist, toponly, includematlab)

% MYDEPFUN creates a dependency matrix
%
% Use as
%    [outlist, depmat] = mydepfun(inlist)
% where
%   inlist   = Nx1 cell-array, descibes the rows
%   outlist  = 1xM cell-array, describes the columns
%   depmat   = NxM matrix, see below
%
% The depmat contains the following values:
%  - 0 if there is no dependency
%  - 1 for the function itself
%  - 2 for a direct dependency

if nargin<2
  toponly = true;
end

if nargin<3
  includematlab = false;
end

for i=1:length(inlist)
  % we first need only the function name
  [p, f, x] = fileparts(inlist{i});
  inlist{i} = f;
  try
    fprintf('processing "%s"\n', inlist{i});
    if toponly
      dep{i} = depfun(inlist{i}, '-toponly', '-quiet'); % determine the dependencies
    else
      dep{i} = depfun(inlist{i},             '-quiet'); % determine the dependencies
    end
    dep{i} = dep{i}(2:end); % the first non-interesting dependency is the function itself
  catch
    dep{i} = {};  % this function cannot be found, so no dependencies
  end
  num(i) = length(dep{i});
  % replace with the full filename, including path
  fullname = which(inlist{i});
  if ~isempty(fullname)
    inlist{i} = fullname;
  else
    inlist{i} = fullfile(p, [f x]);
  end
end

% first make a matrix that will hold all dependencies, including many double ones
depmat = zeros(length(inlist), sum(num));

% mark the dependencies with a value of 2
for i=1:length(inlist)
  offset = sum(num(1:(i-1)));
  depmat(i, (offset+1):(offset+num(i))) = 2;
end
dep(cellfun('isempty', dep))=[];
outlist = cat(1, dep{:});

% add the functions from the inlist as dependent on them selves, i.e. value of 1
depmat  = cat(2, eye(length(inlist)), depmat);
outlist = cat(1, inlist(:), outlist(:));

% remove all double occurences
s = unique(outlist);
for i=1:length(s)
  sel = find(strcmp(outlist, s{i}));
  if numel(sel)>1
    depmat(:,sel(1)) = max(depmat(:,sel),[], 2);
    outlist(sel(2:end))  = {''};  % flag for removal
    depmat(:,sel(2:end)) = nan;   % flag for removal
  end
end
% remove the flagged ones
outlist = outlist(~strcmp(outlist, ''));
depmat  = depmat(:,~any(isnan(depmat),1));

if ~includematlab
  % remove the dependencies on matlab
  sel = false(size(outlist));
  for i=1:length(outlist)
    sel(i) = ~isempty(strfind(outlist{i}, '/opt/matlab')) | ~isempty(strfind(outlist{i}, '/Application'));
  end
  depmat  = depmat(:,~sel);
  outlist = outlist(~sel);
end
