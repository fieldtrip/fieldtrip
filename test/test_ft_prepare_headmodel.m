function test_ft_prepare_headmodel

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_prepare_headmodel

nchan = 64;

% construct a set of electrodes randomly distributed over the upper hemisphere
elec.label = cellstr(num2str((1:nchan).'));
elec.unit = 'cm';
elec.tra = eye(nchan);
elec.elecpos = randn(nchan,3);
elec.elecpos(:,3) = abs(elec.elecpos(:,3));
for i=1:nchan
  elec.elecpos(i,:) = 10*elec.elecpos(i,:)/norm(elec.elecpos(i,:));
end
elec.chanpos = elec.elecpos;

cfg = [];
cfg.method = 'concentricspheres';
cfg.conductivity = [0.33 1.00 0.042 0.33];
headmodel = ft_prepare_headmodel(cfg,elec);
