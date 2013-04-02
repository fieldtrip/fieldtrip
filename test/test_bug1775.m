% function test_bug1775

% TEST test_bug1775
% TEST ft_sourceparcellate ft_datatype_source ft_datatype_parcellation ft_datatype_segmentation



%% create a set of sensors

[pnt, tri] = icosahedron162;

pnt = pnt .* 10; % convert to cm
sel = find(pnt(:,3)>0);

grad.pnt = pnt(sel,:) .* 1.2;
grad.ori = pnt(sel,:);
grad.tra = eye(length(sel));
for i=1:length(sel)
  grad.ori(i,:) = grad.ori(i,:) ./ norm(grad.ori(i,:));
  grad.label{i} = sprintf('magnetometer%d', i);
end
grad.unit = 'cm';
grad.type = 'magnetometer';

grad = ft_datatype_sens(grad);


%% create a volume conductor

vol = [];
vol.r = 10;
vol.o = [0 0 0];
vol.unit = 'cm';

vol = ft_datatype_headmodel(vol);

%% create some precomputed leadfields

cfg = [];
cfg.grad = grad;
cfg.vol = vol;
cfg.resolution = 4;
cfg.channel = 'all';
grid = ft_prepare_leadfield(cfg);

%% create an anatomical parcellation




%% parcellate the leadfields



%% compute a source reconstruction






