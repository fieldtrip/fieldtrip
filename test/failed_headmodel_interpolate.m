function failed_headmodel_interpolate

% WALLTIME 00:45:00
% MEM 6gb

% TEST icosahedron162 ft_voltype ft_headmodel_interpolate ft_prepare_vol_sens ft_compute_leadfield leadfield_interpolate ft_apply_transform

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create a set of electrodes nicely covering the upper half of a sphere
[pnt, tri] = icosahedron162;
pnt = pnt .* 10; % convert to cm
sel = find(pnt(:,3)>0);
elec1 = [];
elec1.pnt = pnt(sel,:);
for i=1:length(sel)
  elec1.label{i} = sprintf('elec%d', i);
end
elec1.unit = 'cm';

% create another set of electrodes covering the upper half of a sphere
% pnt = randn(40,3);
% pnt(:,3) = abs(pnt(:,3)); % only positive z-values
pnt = elec1.pnt(1:42,:) + randn(42,3);
elec2 = [];
for i=1:size(pnt,1)
  elec2.pnt(i,:) = pnt(i,:)./norm(pnt(i,:)) * 10;
  elec2.label{i} = sprintf('rand%d', i);
end
elec2.unit = 'cm';

% update the electrode sets to the latest standards, i.e. elecpos+chanpos rather than pnt
elec1 = ft_datatype_sens(elec1);
elec2 = ft_datatype_sens(elec2);

% create another set of electrodes representing a bipolar montage
bipolar = [];
bipolar.labelold = elec1.label;
bipolar.tra = zeros(length(elec1.label)-1, length(elec1.label));
for i=1:(length(bipolar.labelold)-1)
  bipolar.labelnew{i} = sprintf('%s-%s', bipolar.labelold{i}, bipolar.labelold{i+1});
  bipolar.tra(i,i  ) = +1;
  bipolar.tra(i,i+1) = -1;
end
elec3 = ft_apply_montage(elec1, bipolar);

% create an identical set of electrodes but with other channel names
elec4 = elec1;
for i=1:length(elec4.label)
  elec4.label{i} = sprintf('chan%d', i);
end

% create another bipolar set of channels that is incompatible with the previous one
bipolar = [];
bipolar.labelold = elec1.label;
bipolar.tra      = zeros(10, length(elec1.label));
for i=1:(10)
  bipolar.labelnew{i} = sprintf('%s-%s', bipolar.labelold{i}, bipolar.labelold{i+2});
  bipolar.tra(i,i  ) = +1;
  bipolar.tra(i,i+2) = -1;
end
elec5 = ft_apply_montage(elec1, bipolar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% construct a singlesphere volume conduction model
volA   = [];
volA.c = 1;
volA.r = 10;
volA.o = [0 0 0];

cfg      = [];
cfg.vol = volA;
cfg.elec = elec1;
cfg.grid.resolution = 1;
leadfield = ft_prepare_leadfield(cfg);

% remember one position
pos1 = leadfield.pos(leadfield.inside(1),:);
pos2 = [0 0 7];

lf1  = leadfield.leadfield{leadfield.inside(1)};

% the following call generates the volume conduction model and writes all files to disk
filename = fullfile(tempname, 'leadfield');
ft_headmodel_interpolate(filename, elec1, leadfield, 'smooth', false);

% the next day you would start by reading it from disk
volB = ft_read_vol([filename '.mat']); % this is a mat file containing the "vol" structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use the same electrodes
[volAA, elecAA] = ft_prepare_vol_sens(volA, elec1);
[volBB, elecBB] = ft_prepare_vol_sens(volB, elec1);

% compare the leadfields
lfa = ft_compute_leadfield(pos1, elecAA, volAA); % the original
lfb = ft_compute_leadfield(pos1, elecBB, volBB); % the interpolation

assert(isalmostequal(lf1, lfb, 'reltol', 1e-4), 'the leadfields are different');
assert(isalmostequal(lfa, lfb, 'reltol', 1e-4), 'the leadfields are different');
assert(all(mean(lfa, 1)<eps(double(1))), 'the leadfield is not average referenced'); % this is in double precision
assert(all(mean(lfb, 1)<eps(single(1))), 'the leadfield is not average referenced'); % this is in single precision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use a subset of electrodes
[volAA, elecAA] = ft_prepare_vol_sens(volA, elec1, 'channel', elec1.label(1:50));
[volBB, elecBB] = ft_prepare_vol_sens(volB, elec1, 'channel', elec1.label(1:50));

% compare the leadfields
lfa = ft_compute_leadfield(pos1, elecAA, volAA); % the original
lfb = ft_compute_leadfield(pos1, elecBB, volBB); % the interpolation

assert(isalmostequal(lfa, lfb, 'reltol', 1e-3), 'the leadfields are different');
assert(all(mean(lfa, 1)<eps), 'the leadfield is not average referenced');
assert(all(mean(lfb, 1)<eps), 'the leadfield is not average referenced');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use the same electrodes with with different names
[volAA, elecAA] = ft_prepare_vol_sens(volA, elec4);
[volBB, elecBB] = ft_prepare_vol_sens(volB, elec4);

% compare the leadfields
lfa = ft_compute_leadfield(pos1, elecAA, volAA); % the original
lfb = ft_compute_leadfield(pos1, elecBB, volBB); % the interpolation

assert(isalmostequal(lfa, lfb, 'reltol', 1e-4), 'the leadfields are different');
assert(all(mean(lfa, 1)<eps), 'the leadfield is not average referenced');
assert(all(mean(lfb, 1)<eps), 'the leadfield is not average referenced');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use electrodes with a different placement
[volAA, elecAA] = ft_prepare_vol_sens(volA, elec2);
[volBB, elecBB] = ft_prepare_vol_sens(volB, elec2);

% compare the leadfields
lfa = ft_compute_leadfield(pos2, elecAA, volAA); % the original
lfb = ft_compute_leadfield(pos2, elecBB, volBB); % the interpolation

c = corrcoef(lfa(:), lfb(:)); c = c(1,2);
assert(c>0.8, 'the leadfields are different');
assert(all(mean(lfa, 1)<eps), 'the leadfield is not average referenced');
assert(all(mean(lfb, 1)<eps), 'the leadfield is not average referenced');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use a bipolar electrode montage
[volAA, elecAA] = ft_prepare_vol_sens(volA, elec3);
[volBB, elecBB] = ft_prepare_vol_sens(volB, elec3);

% compare the leadfields
lfa = ft_compute_leadfield(pos2, elecAA, volAA); % the original
lfb = ft_compute_leadfield(pos2, elecBB, volBB); % the interpolation

assert(isalmostequal(lfa, lfb, 'abstol', inf, 'reltol', inf, 'relnormtol', 1e-3), 'the leadfields are different');
assert(~all(mean(lfa, 1)<eps), 'the leadfield is average referenced');
assert(~all(mean(lfb, 1)<eps), 'the leadfield is average referenced');
