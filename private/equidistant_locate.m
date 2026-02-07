function [elc, lab] = equidistant_locate(pos, tri, Fpz, Oz, T7, T8, numelec, nummidline, feedback)

% EQUIDISTANT_LOCATE - Place electrodes equidistantly on the mesh surface
%
% See also ELEC1020_LOCATE

if nargin<8
  feedback = true;
end

if feedback
  figure
end

% to avoid confusion between electrode and headshape positions
headshape.pos = pos;
headshape.tri = tri;
clear pos tri
headshape.nrm = surface_normals(headshape.pos, headshape.tri);

% estimate how large the head is
headsize = (norm(T7-T8) + norm(Fpz - Oz)) / 2;

% determine the midpoint of the fiducials and the cardinal directions
midpoint = (T7 + T8)/2;
anterior = Fpz - midpoint;
anterior = anterior / norm(anterior);
superior = cross(T8 - T7, anterior);
superior = superior / norm(superior);
right = cross(anterior, superior);
right = right / norm(right);

Fpz = Fpz - midpoint;
Oz  = Oz  - midpoint;
T7  = T7  - midpoint;
T8  = T8  - midpoint;

headshape.pos(:,1) = headshape.pos(:,1) - midpoint(1);
headshape.pos(:,2) = headshape.pos(:,2) - midpoint(2);
headshape.pos(:,3) = headshape.pos(:,3) - midpoint(3);

% create some electrodes along the midline
midlineelc = zeros(nummidline, 3);


% crete some electrodes on the left hemisphere
leftelec = zeros((numelec-nummidline-4)/2, 3);


elc = randn(numelec-4,3);
elc(2:2:end,1) = -elc(1:2:end,1);
elc(2:2:end,2) =  elc(1:2:end,2);
elc(:,3) = sqrt(2)/2;
if mod(numelec-4,2)
  elc(end,1) = 0;
  elc(end,2) = 0;
end
for i=1:size(elc,1)
  elc(i,:) = elc(i,:) / norm(elc(i,:));
end
elc = elc * headsize/2;


% add the fixed electrodes
elc = [
  elc
  Fpz
  Oz
  T7
  T8
  ];

% project the electrodes onto the mesh surface
[dum, elc] = project_elec(elc, headshape.pos, headshape.tri);

stepsize = headsize/1000;
minchange = headsize/10000;
maxiter   = 1000;

change = inf;
iter = 0;
while change>minchange && iter<maxiter
  % take the previous positions
  pos = elc;

  % make a triangulation of the electrodes
  prj = elproj(pos);
  tri = delaunay(prj(:,1), prj(:,2));

  edge = [
    tri(:,[1 2])
    tri(:,[2 3])
    tri(:,[3 1])
    ];

  % compute the length of each edge
  l = sqrt(sum((pos(edge(:,1),:)-pos(edge(:,2),:)).^2,2));
  s = std(l);

  % loop over the edges and for each edge compute how much the two vertices that it
  % connects should move towards or away from each other 
  delta = zeros(size(pos));
  for i=1:size(edge,1)
    v1 = edge(i,1);
    v2 = edge(i,2);
    ori = pos(v1,:) - pos(v2,:);
    ori = ori / norm(ori);

    sel = any(edge==v1, 2); % selection of the neighbours
    m = mean(l(sel));       % mean distance to the neighbours
    k = sum(sel);           % number of neighbours
    delta(v1,:) = delta(v1,:) + (stepsize/k) * (m - l(i)) * ori;

    sel = any(edge==v1, 2); % selection of the neighbours
    m = mean(l(sel));       % mean distance to the neighbours
    k = sum(sel);           % number of neighbours
    delta(v2,:) = delta(v2,:) + (stepsize/k) * (l(i) - m) * ori;
  end

  for i=1:numelec
    if delta(i,3)<0
      % they are not allowed to shift below the fixed electrodes
      delta(i,3) = 0;
    end
  end

  % the last four are the reference electrodes which are not to move in any case
  delta(end-3,:) = 0;
  delta(end-2,:) = 0;
  delta(end-1,:) = 0;
  delta(end-0,:) = 0;

  % update the positions
  pos = pos + delta;

  % project the new positions onto the mesh surface
  [dum, pos] = project_elec(pos, headshape.pos, headshape.tri);

  % compute the average shift of the electrodes
  change = mean(sqrt(sum((elc - pos).^2,2)));
  iter = iter + 1;

  if feedback
    cla
    ft_plot_mesh(struct('pos', pos, 'tri', tri))
    drawnow
    fprintf('iter = %d, mean = %f, std = %f, change = %f\n', iter, m, s, change)
  end

  % repeat with the updated positions
  elc = pos;
end % while

lab = cell(numelec, 1);
for i=1:numelec
  lab{i} = sprintf('%d', i);
end
% the last four electrodes are the fixed ones
lab{n-3} = 'Fpz';
lab{n-2} = 'Oz';
lab{n-1} = 'T7';
lab{n-0} = 'T8';
