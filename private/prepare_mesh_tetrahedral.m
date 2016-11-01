function mesh = prepare_mesh_tetrahedral(cfg, mri)

% ensure that the input is consistent with what this function expects
mri = ft_checkdata(mri, 'datatype', {'volume', 'segmentation'}, 'hasunit', 'yes');

% get the default options
cfg.tissue = ft_getopt(cfg, 'tissue');

if isempty(cfg.tissue)
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
  fn = fieldnames(mri);
  for i = 1:numel(fn)
    if (numel(mri.(fn{i})) == prod(mri.dim)) && (~strcmp(fn{i}, 'inside'))
      segfield = fn{i};
    end
  end
  cfg.tissue = setdiff(unique(mri.(segfield)(:)), 0);
end

if ischar(cfg.tissue)
  % it should either be something like {'brain', 'skull', 'scalp'}, or something like [1 2 3]
  cfg.tissue = {cfg.tissue};
end

if iscell(cfg.tissue)
  % the code below assumes that it is a probabilistic representation
  if any(strcmp(cfg.tissue, 'brain'))
    mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic', 'hasbrain', 'yes');
  else
    mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic');
  end
  % combine all tissue types
  seg = false(mri.dim);
  for i=1:numel(cfg.tissue)
    seg = seg | mri.(cfg.tissue{i});
  end
else
  % the code below assumes that it is an indexed representation
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
  % combine all tissue types
  seg = (mri.seg>0);
end

isovalue = 0.5;

[node, elem, face] = vol2mesh(seg, 1:mri.dim(1), 1:mri.dim(2), 1:mri.dim(3), 2, 2, isovalue);

mesh = keepfields(mri, {'coordsys', 'unit'});
mesh.pos = ft_warp_apply(mri.transform, node);
mesh.tet = elem(:,1:4);