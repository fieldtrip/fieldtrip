function test_ft_prepare_leadfield

% WALLTIME 00:20:00
% MEM 2gb
% DEPENDENCY ft_prepare_leadfield

fs = 500;
nchan = 32;
start_time = -1; % seconds
end_time = 2.5; % seconds
nsamples = (end_time - start_time) * fs + 1;

data = [];
data.time{1} = linspace(start_time, end_time, nsamples);
data.trial{1} = randn(nchan,nsamples);
data.sampleinfo = [1 nsamples];
data.label = cellstr(num2str((1:nchan).'));

% construct a set of electrodes randomly distributed over the upper hemisphere
data.elec.label = data.label;
data.elec.unit = 'cm';
data.elec.tra = eye(nchan);
data.elec.elecpos = randn(nchan,3);
data.elec.elecpos(:,3) = abs(data.elec.elecpos(:,3));
for i=1:nchan
  data.elec.elecpos(i,:) = 10*data.elec.elecpos(i,:)/norm(data.elec.elecpos(i,:));
end
data.elec.chanpos = data.elec.elecpos;

% make an explicit geometry with 4 surfaces 
% the radius of each follows the BESA defaults
geometry = [];
geometry(1).pos  = data.elec.elecpos*71/85; % brain
geometry(1).unit = data.elec.unit;
geometry(2).pos  = data.elec.elecpos*72/85; % csf
geometry(2).unit = data.elec.unit;
geometry(3).pos  = data.elec.elecpos*79/85; % skull
geometry(3).unit = data.elec.unit;
geometry(4).pos  = data.elec.elecpos*85/85; % scalp
geometry(4).unit = data.elec.unit;

% fit a 4-sphere concentric model to the geometry
headmodel = ft_headmodel_concentricspheres(geometry, 'conductivity', [0.33 1.00 0.042 0.33]);

cfg = [];
cfg.channel = data.label;
cfg.elec = data.elec;
cfg.headmodel = headmodel;
sourcemodel = ft_prepare_leadfield(cfg, data);
