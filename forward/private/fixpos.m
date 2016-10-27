function mesh = fixpos(mesh, recurse)

% FIXPOS helper function to ensure that meshes are described properly

if nargin==1
  recurse = 1;
end

if isa(mesh, 'delaunayTriangulation')
  % convert to structure, otherwise the code below won't work properly
  ws = warning('off', 'MATLAB:structOnObject');
  mesh = struct(mesh);
  warning(ws);
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

% convert from MATLAB delaunayTriangulation output to FieldTrip convention
if isfield(mesh, 'Points') && isfield(mesh, 'ConnectivityList')
  mesh.pos = mesh.Points;
  switch size(mesh.ConnectivityList,2)
    case 2
      mesh.edge = mesh.ConnectivityList;
    case 3
      mesh.tri = mesh.ConnectivityList;
    case 4
      mesh.tet = mesh.ConnectivityList;
    case 8
      mesh.hex = mesh.ConnectivityList;
    otherwise
      error('unsupported ConnectivityList')
  end % switch
  mesh = removefields(mesh, {'Points', 'ConnectivityList', 'Constraints', 'UnderlyingObj'});
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
