function montage = fixmontage(montage, revert)

% FIXMONTAGE use "old/new" instead of "org/new"

if nargin<2
  revert = false;
end

if ~revert
  from = {'labelorg', 'chantypeorg', 'chanunitorg', 'chanposorg', 'chanoriorg'};
  to   = {'labelold', 'chantypeold', 'chanunitold', 'chanposold', 'chanoriold'};
else
  to   = {'labelorg', 'chantypeorg', 'chanunitorg', 'chanposorg', 'chanoriorg'};
  from = {'labelold', 'chantypeold', 'chanunitold', 'chanposold', 'chanoriold'};
end

for i=1:numel(from)
  if isfield(montage, from{i})
    montage.(to{i}) = montage.(from{i});
    montage = rmfield(montage, from{i});
  end
end

% use recursion to update the subfields
fn = fieldnames(montage);
for i=1:numel(fn)
  if isstruct(montage.(fn{i}))
    montage.(fn{i}) = fixmontage(montage.(fn{i}), revert);
  end
end
