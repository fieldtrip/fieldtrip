function test_bug3453

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY

%%

% unit radius, scaling is done further down
pos = randn(100,3);
for i=1:100
  pos(i,:) = pos(i,:)/norm(pos(i,:));
end
% brain, csf, skull, skin
headshape = [];
headshape.bnd(1).unit = 'm';
headshape.bnd(2).unit = 'm';
headshape.bnd(3).unit = 'm';
headshape.bnd(4).unit = 'm';
headshape.bnd(1).pos = pos*0.087;
headshape.bnd(2).pos = pos*0.088;
headshape.bnd(3).pos = pos*0.092;
headshape.bnd(4).pos = pos*0.100;

elec = [];
elec.unit = 'm';
elec.elecpos = randn(40,3);
for i=1:40
  elec.label{i} = num2str(i);
  elec.elecpos(i,:) = 0.1*elec.elecpos(i,:)/norm(elec.elecpos(i,:));
end
elec.chanpos = elec.elecpos;
elec.tra = eye(40);

dippos = [0 0 0.08];

%%

cfg = [];
cfg.method = 'concentricspheres';
cfg.headshape = headshape;
cfg.conductivity = [0.3300 1 0.0042 0.3300]; % brain, csf, skull, skin
vol1 = ft_prepare_headmodel(cfg);

cfg.order = 60;
vol2 = ft_prepare_headmodel(cfg); % should be same as default

cfg.order = 5;
vol3 = ft_prepare_headmodel(cfg); % very inaccurate

%% using high-level code

cfg = [];
cfg.sourcemodel.pos = dippos;
cfg.sourcemodel.unit = 'm';
cfg.elec = elec;
cfg.headmodel = vol1;
lf1 = ft_prepare_leadfield(cfg);
cfg.headmodel = vol2;
lf2 = ft_prepare_leadfield(cfg);
cfg.headmodel = vol3;
lf3 = ft_prepare_leadfield(cfg);

%% slightly lower level code

% note that I am skipping ft_prepare_vol_sens, which would not make a difference here
lf1 = ft_compute_leadfield(dippos, elec, vol1);
lf2 = ft_compute_leadfield(dippos, elec, vol2);
lf3 = ft_compute_leadfield(dippos, elec, vol3); % should be inaccurate

assert( isequal(lf1, lf2));
assert(~isequal(lf1, lf3));
