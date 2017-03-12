function input = fixoldorg(input, revert)

% FIXOLDORG use "old/new" instead of "org/new"

if nargin<2
  revert = false;
end

if ~(isstruct(input) && numel(input)==1)
  % this does not look like a montage or a sensor description
  return
end

if ~revert
  from = {'labelorg', 'chantypeorg', 'chanunitorg', 'chanposorg', 'chanoriorg'};
  to   = {'labelold', 'chantypeold', 'chanunitold', 'chanposold', 'chanoriold'};
else
  to   = {'labelorg', 'chantypeorg', 'chanunitorg', 'chanposorg', 'chanoriorg'};
  from = {'labelold', 'chantypeold', 'chanunitold', 'chanposold', 'chanoriold'};
end

% replace the fields
for i=1:numel(from)
  if isfield(input, from{i})
    input.(to{i}) = input.(from{i});
    input = rmfield(input, from{i});
  end
end

% use recursion to update the subfields
fn = fieldnames(input);
for i=1:numel(fn)
  if isstruct(input.(fn{i}))
    input.(fn{i}) = fixoldorg(input.(fn{i}), revert);
  end
end
