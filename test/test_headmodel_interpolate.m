% function test_headmodel_interpolate

% TEST test_headmodel_interpolate
% TEST icosahedron162 ft_voltype ft_headmodel_interpolate ft_prepare_vol_sens ft_compute_leadfield leadfield_interpolate

[pnt, tri] = icosahedron162;

% create sensors in cm
pnt = pnt .* 10; % convert to cm
sel = find(pnt(:,3)>0);

elec.pnt = pnt(sel,:);
for i=1:length(sel)
  elec.label{i} = sprintf('electrode%d', i);
end
elec.unit = 'cm';

othersens = [];
pnt = randn(40,3);
pnt(:,3) = abs(pnt(:,3)); % only positive z-values
for i=1:40
  pnt(i,:) = pnt(i,:)./norm(pnt(i,:));
end

% update it to the latest standards
elec = ft_datatype_sens(elec);

vol1   = [];
vol1.c = 1;
vol1.r = 10;
vol1.o = [0 0 0];

cfg      = [];
cfg.vol  = vol1;
cfg.elec = elec;
cfg.grid.resolution = 2;
leadfield = ft_prepare_leadfield(cfg);

filename = tempname;

ft_headmodel_interpolate(filename, elec, leadfield);

% the next day you would do
vol = ft_read_vol(filename); % this is a mat file containing a "vol" structure

[vol, sens] = ft_prepare_vol_sens(vol, sens);
lf = ft_compute_leadfield(randn(10,3), sens, vol);

[vol, sens] = ft_prepare_vol_sens(vol, sens, 'channel', sens.label(1:10));
lf = ft_compute_leadfield(randn(10,3), sens, vol);

[vol, othersens] = ft_prepare_vol_sens(vol, othersens);
lf = ft_compute_leadfield(randn(10,3), sens, vol);






