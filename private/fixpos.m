function data = fixpos(data)

% helper function that replaces pnt by pos

if isfield(data, 'pnt')
  data.pos = data.pnt;
  data = rmfield(data, 'pnt');
end

% recurse into substructures
fn = fieldnames(data);
fn = setdiff(fn, {'cfg'});
for i=1:length(fn)
  if isstruct(data.(fn{i}))
    data.(fn{i}) = fixpos(data.(fn{i}));
  end
end