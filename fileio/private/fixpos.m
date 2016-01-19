function mesh = fixpos(mesh, recurse)

% helper function to replace pnt by pos

if nargin==1
  recurse = 1;
end

if ~isa(mesh, 'struct')
    return;
end

if numel(mesh)>1
  % loop over all individual elements
  clear tmp
  for i=1:numel(mesh)
    % this is to prevent an "Subscripted assignment between dissimilar structures" error
    tmp(i) = fixpos(mesh(i));
  end
  mesh = tmp;
  clear tmp
  return
end

% convert from BrainStorm/MNE to FieldTrip convention
if isfield(mesh, 'vertices') && isfield(mesh, 'faces')
  mesh.pos = mesh.vertices;
  mesh.tri = mesh.faces;
  mesh = rmfield(mesh, {'faces', 'vertices'});
elseif isfield(mesh, 'Vertices') && isfield(mesh, 'Faces')
  mesh.pos = mesh.Vertices;
  mesh.tri = mesh.Faces;
  mesh = rmfield(mesh, {'Faces', 'Vertices'});
end
  
% replace pnt by pos
if isfield(mesh, 'pnt')
  mesh.pos = mesh.pnt;
  mesh = rmfield(mesh, 'pnt');
end

if recurse<3
  % recurse into substructures, not too deep
  fn = fieldnames(mesh);
  fn = setdiff(fn, {'cfg'}); % don't recurse into the cfg structure
  for i=1:length(fn)
    if isstruct(mesh.(fn{i}))
      mesh.(fn{i}) = fixpos(mesh.(fn{i}), recurse+1);
    end
  end
end
