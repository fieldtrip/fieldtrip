function test_ft_dipolefitting_hartmut

% WALLTIME 00:30:00
% MEM 2gb
% DEPENDENCY ft_dipolefitting ft_prepare_sourcemodel ft_inside_headmodel ft_inverse_dipolefit hartmut_eyemodel hartmut_add_eyes
% DATA no

% This tests the HArtMuT extension of ft_dipolefitting: fitting sources in the scalp
% compartment (muscle artefacts) and fitting the eyes as a mirror-symmetric dipole pair.
% See ft_dipolefitting (cfg.dipfit.hartmut) and the design notes for the rationale.
%
% Checks 1 to 4 are self-contained and need no external solver. Checks 5 to 7 need real
% EEG leadfields from a BEM, so they are gated behind OpenMEEG and skipped when it is
% missing. Nothing here depends on private data files; an optional pin against the shipped
% extended_sourcemodel.mat is documented in check 4 but not required to pass.

rng(42); % deterministic synthetic data

hasopenmeeg = ft_hastoolbox('openmeeg', 3); % 3 = check silently, do not error
ftpath = fileparts(which('ft_dipolefitting')); % to reach private functions, see check 4

%% check 1: regression, hartmut off leaves a brain fit unchanged
% A brain dipole fit must give the same result whether cfg.dipfit.hartmut is absent or 'no'.
% This proves the new config path is inert when the extension is off.

[elec, sphere] = local_spheremodel();

% simulate a single brain dipole and average it to a timelock structure
truepos = [2 3 6]; % cm, inside the brain sphere
cfgs = [];
cfgs.headmodel = sphere;
cfgs.elec      = elec;
cfgs.dip.pos   = truepos;
cfgs.dip.mom   = [1 0 0]';
cfgs.dip.frequency = 10;
cfgs.dip.phase     = 0;
cfgs.fsample   = 250;
cfgs.relnoise  = 0;
cfgs.ntrials   = 1;
sim = ft_dipolesimulation(cfgs);
timelock = ft_timelockanalysis([], sim);

cfgf = [];
cfgf.numdipoles = 1;
cfgf.headmodel  = sphere;
cfgf.elec       = elec;
cfgf.gridsearch = 'yes';
cfgf.nonlinear  = 'yes';
cfgf.resolution = 2; % cm, coarse for speed
cfgf.model      = 'regional';
cfgf.latency    = 'all';

% reference fit, no hartmut field at all
src_ref = ft_dipolefitting(cfgf, timelock);

% same fit with the extension explicitly off
cfgf.dipfit.hartmut     = 'no';
cfgf.dipfit.compartment = 'brain';
src_off = ft_dipolefitting(cfgf, timelock);

assert(isequal(src_ref.dip.pos, src_off.dip.pos), 'hartmut=no changed the fitted position');
assert(isequal(src_ref.dip.mom, src_off.dip.mom), 'hartmut=no changed the fitted moment');
% the fit should also be near the simulated source
assert(norm(src_ref.dip.pos - truepos) < 2, 'the brain dipole was not recovered'); % cm

%% check 2: ft_inside_headmodel labels the scalp compartment correctly
% A point between skin and skull is in the scalp; points in skull, csf, brain, and outside
% the skin are not. This is pure geometry and needs no solver.

bem = local_nestedspheres(); % bnd ordered outside-in, type openmeeg, mm

rscalp = 90; rskull = 82; % see local_nestedspheres
inscalp  = [0 0 86]; % between skin (90) and skull (82)
inskull  = [0 0 80]; % inside the skull boundary
inbrain  = [0 0 10]; % deep inside
outside  = [0 0 120]; % beyond the skin

testpos = [inscalp; inskull; inbrain; outside];
inside  = ft_inside_headmodel(testpos, bem, 'surface', 'scalp');
assert(inside(1)==true,  'a point between skin and skull should be in the scalp');
assert(inside(2)==false, 'a point inside the skull should not be in the scalp');
assert(inside(3)==false, 'a deep brain point should not be in the scalp');
assert(inside(4)==false, 'a point outside the skin should not be in the scalp');

% the brain compartment still works as before
insidebrain = ft_inside_headmodel(testpos, bem, 'surface', 'brain');
assert(insidebrain(3)==true,  'a deep point should be in the brain');
assert(insidebrain(1)==false, 'a scalp point should not be in the brain');

% a model with a single boundary cannot define a scalp shell
onebnd = bem; onebnd.bnd = bem.bnd(1);
try
  ft_inside_headmodel(inscalp, onebnd, 'surface', 'scalp');
  error('expected an error for a scalp test on a single-boundary model');
catch me
  assert(~isempty(strfind(me.message, 'two boundaries')) || ~isempty(strfind(lower(me.message), 'boundaries')), ...
    'wrong error for a single-boundary scalp test');
end

%% check 3: linkmom recovers a shared moment for a symmetric pair
% Simulate data from one moment shared by a mirrored pair (mom2 = mommap*mom1), fit with
% the linkmom constraint, and confirm the moment comes back and stays linked.

mommap = eye(3); % M = I, synchronous eyes
p1 = [3 5 4]; % cm, inside the brain so checkinside is irrelevant
p2 = p1 .* [-1 1 1]; % mirror across the midsagittal (x) plane

lf1 = ft_compute_leadfield(p1, elec, sphere);
lf2 = ft_compute_leadfield(p2, elec, sphere);
m1  = [0.7 -0.3 0.5]'; % the true shared moment
Vdata = lf1*m1 + lf2*(mommap*m1); % equals (lf1 + lf2*M)*m1

dipin = [];
dipin.pos = [p1; p2];
dipin.mom = zeros(6,1);

constr = [];
constr.reduce  = [1 2 3];
constr.expand  = [1 2 3 1 2 3];
constr.mirror  = ones(1,6); constr.mirror(3+1) = -1; % flip x of the second dipole (setting mirror sets symmetry=true)
constr.linkmom = true;
constr.mommap  = mommap;

est = ft_inverse_dipolefit(dipin, elec, sphere, Vdata, 'constr', constr, 'checkinside', false);

mom = est.mom(:);
assert(numel(mom)==6, 'a linkmom fit must return a 6x1 moment for the pair');
% the link itself is exact: the second moment is always mommap times the first
assert(norm(mom(4:6) - mommap*mom(1:3)) < 1e-6, 'the second moment is not mommap*first');
% the two dipoles stay exact mirror images, which confirms the symmetric reduction
assert(norm(est.pos(2,:) - est.pos(1,:).*[-1 1 1]) < 1e-6, 'the pair is not symmetric across x');
% the recovered moment and position match the simulated source up to the EEG average
% reference convention. ft_inverse_dipolefit average references the data but builds the
% model from a raw leadfield, so the true source is not the exact least squares minimum and
% a few percent of moment error and a few mm of position drift are expected, the same as for
% an ordinary unconstrained single dipole fit. The exact checks above pin the linkmom wiring.
assert(norm(mom(1:3) - m1)/norm(m1) < 0.1, 'the shared moment was not recovered');
assert(norm(est.pos(1,:) - p1) < 0.5, 'the symmetric pair did not land near the true source'); % cm

% the linkmom constraint must reject incompatible options
try
  badc = constr; badc.fixedori = true;
  ft_inverse_dipolefit(dipin, elec, sphere, Vdata, 'constr', badc, 'checkinside', false);
  error('expected linkmom to reject fixedori');
catch me
  assert(~isempty(strfind(lower(me.message), 'linkmom')), 'wrong error for linkmom with fixedori');
end

%% check 4: the fused eye leadfield follows the M=I convention
% hartmut_add_eyes stores each ocular candidate as L(p) + L(mirror(p))*mommap. Confirm that
% the stored fused leadfield equals L(p) + L(mirror(p)) with mommap = I, which pins the
% convention used everywhere downstream.
%
% Note on the shipped file: raw/dipfit/extended_BEM/extended_sourcemodel.mat carries 12
% symmetric_dipole eye entries with the same convention. Reproducing them needs the exact
% BEM and OpenMEEG that built the file, so that pin is an offline check, not part of CI. The
% convention itself is what we test here, against the same code path the file was built with.

eyesphere = sphere;
eyesphere.coordsys = 'mni'; % so hartmut_eyemodel uses its default MNI eye geometry
eyesphere.unit     = 'cm';
elecmni = elec; elecmni.coordsys = 'mni';

% a near-empty source model; hartmut_add_eyes only appends the eye candidates
sm = [];
sm.pos    = zeros(0,3);
sm.unit   = 'cm';
sm.inside = false(0,1);

cfge = [];
cfge.dipfit.eye.radius      = 8;  % mm, small so the candidate grid stays tiny
cfge.dipfit.eye.interocular = 70; % mm
cfge.dipfit.eye.offset      = [72 -25]; % mm
cfge.dipfit.constr.mommap   = eye(3);

leadfieldopt = {}; % defaults

% hartmut_add_eyes lives in fieldtrip/private, so cd into it the way other FieldTrip tests
% reach private functions, then restore the working directory
here = pwd;
cleanup = onCleanup(@() cd(here));
cd(fullfile(ftpath, 'private'));
sm = hartmut_add_eyes(sm, cfge, elecmni, eyesphere, leadfieldopt);
clear cleanup; % triggers the cd back

assert(~isempty(sm.pos), 'hartmut_add_eyes added no ocular candidates');
assert(all(strcmp(sm.compartment, 'eye')), 'the added candidates are not labelled as eye');

% recompute the fused leadfield independently and compare
mirror = [-1 1 1]; % x is the left-right axis in MNI
for i=1:size(sm.pos,1)
  p   = sm.pos(i,:);
  lfa = ft_compute_leadfield(p,          elecmni, eyesphere);
  lfb = ft_compute_leadfield(p.*mirror,  elecmni, eyesphere);
  expected = lfa + lfb*eye(3);
  assert(norm(sm.leadfield{i} - expected, 'fro') < 1e-9, 'the stored eye leadfield breaks the M=I fuse convention');
end

%% checks 5 to 7 need OpenMEEG for real EEG leadfields
if ~hasopenmeeg
  ft_warning('OpenMEEG is not available, skipping the scalp, eye, and end-to-end fits (checks 5 to 7)');
  return
end

% build an OpenMEEG BEM from the nested spheres, in mm and MNI-like for the eye geometry
ommodel = ft_headmodel_openmeeg(bem.bnd, 'conductivity', [0.43 0.01 1.79 0.33]);
ommodel.coordsys = 'mni';
ommodel.unit     = 'mm';
elecmm = local_projectelec(elec, 90); % put electrodes on the 90 mm scalp surface
elecmm.coordsys = 'mni';

%% check 5: a muscle dipole in the scalp is recovered in the scalp
% fit one muscle topography as a single component, the same way a real ICA component is fit
% (see checks 6 and 7).
scalppos = [0 50 70]; % mm, norm 86, between the skull (82) and the skin (90)
lfScalp  = ft_compute_leadfield(scalppos, elecmm, ommodel);
compS    = local_topo2comp(lfScalp*[0 0 1]', elecmm);

cfgf = [];
cfgf.headmodel  = ommodel;
cfgf.elec       = elecmm;
cfgf.gridsearch = 'yes';
cfgf.nonlinear  = 'yes';
cfgf.resolution = 10; % mm, fine enough to seed the thin scalp shell
cfgf.model      = 'moving';
cfgf.component  = 1;
cfgf.dipfit.hartmut     = 'yes';
cfgf.dipfit.compartment = 'scalp'; % look in the scalp only for this check
srcS = ft_dipolefitting(cfgf, compS);

fitpos = srcS.dip(1).pos(1,:);
assert(ft_inside_headmodel(fitpos, ommodel, 'surface', 'scalp'), 'the muscle source was not fit inside the scalp');

%% check 6: a symmetric ocular source comes back as a linked dipole pair
% simulate a synchronous pair at the default MNI eye centres, then fit moving. The default
% centre is [-interocular/2, offset(1), offset(2)] = [-35 72 -25] mm, see hartmut_eyemodel;
% hartmut_add_eyes (exercised in check 4) builds it the same way during the fit.
pL = [-35 72 -25];          % left eye centre, mm
pR = pL .* [-1 1 1];        % right eye centre
mom = [0.8 0.2 -0.4]';

lfL = ft_compute_leadfield(pL, elecmm, ommodel);
lfR = ft_compute_leadfield(pR, elecmm, ommodel);
Veye = lfL*mom + lfR*mom; % synchronous, M = I

comp = local_topo2comp(Veye, elecmm);

cfgf = [];
cfgf.headmodel  = ommodel;
cfgf.elec       = elecmm;
cfgf.gridsearch = 'yes';
cfgf.nonlinear  = 'yes';
cfgf.resolution = 10; % mm
cfgf.model      = 'moving';
cfgf.component  = 1;
cfgf.dipfit.hartmut = 'yes';
srcE = ft_dipolefitting(cfgf, comp);

eyedip = srcE.dip(1);
assert(size(eyedip.pos,1)==2, 'the ocular component should be a 2-dipole pair');
me = eyedip.mom(:);
assert(norm(me(4:6) - me(1:3)) < 1e-6*max(1,norm(me(1:3))), 'the eye moments are not linked (M=I)');
% the two dipoles must be mirror images across the midsagittal (x) plane
assert(norm(eyedip.pos(2,:) - eyedip.pos(1,:).*[-1 1 1]) < 1e-6, 'the eye pair is not symmetric across x');
% the structure above (a symmetric, linked-moment pair) is the part that is pinned exactly.
% the exact depth of the pair is biased by the EEG average reference convention:
% ft_inverse_dipolefit average references the data while the model leadfield stays raw,
% shifting the optimum by a few cm. localisation is therefore only asserted loosely here.
d2eye = min(norm(eyedip.pos(1,:) - pL), norm(eyedip.pos(1,:) - pR));
assert(d2eye < 45, 'the eye dipole drifted unreasonably far from the eyes'); % mm, loose, see PR notes

%% check 7: end to end, brain + muscle + ocular in one moving fit
brainpos = [20 30 60]; momB = [1 0 0]'; % mm, norm 70, inside the brain (78)
musclepos = [75 0 40]; momM = [0 0 1]'; % mm, norm 85, lateral scalp (temporal muscle)

lfB = ft_compute_leadfield(brainpos, elecmm, ommodel);
lfM = ft_compute_leadfield(musclepos, elecmm, ommodel);
Vb = lfB*momB;
Vm = lfM*momM;
% reuse the ocular topography Veye from check 6
comp3 = local_topo2comp([Vb Vm Veye], elecmm);

cfgf = [];
cfgf.headmodel  = ommodel;
cfgf.elec       = elecmm;
cfgf.gridsearch = 'yes';
cfgf.nonlinear  = 'yes';
cfgf.resolution = 10; % mm
cfgf.model      = 'moving';
cfgf.component  = [1 2 3];
cfgf.dipfit.hartmut = 'yes';
src3 = ft_dipolefitting(cfgf, comp3);

% component 1 brain, component 2 scalp muscle, component 3 symmetric eyes
assert(ft_inside_headmodel(src3.dip(1).pos(1,:), ommodel, 'surface', 'brain'), 'component 1 was not fit in the brain');
assert(ft_inside_headmodel(src3.dip(2).pos(1,:), ommodel, 'surface', 'scalp'), 'component 2 was not fit in the scalp');
assert(size(src3.dip(3).pos,1)==2, 'component 3 should be a symmetric eye pair');

end % main function


%% local helpers

function [elec, headmodel] = local_spheremodel()
% a 4-shell concentric sphere EEG model in cm, with electrodes on the upper hemisphere
nchan = 32;
elec = [];
elec.label = cellstr(num2str((1:nchan)'));
elec.unit  = 'cm';
elec.elecpos = randn(nchan,3);
elec.elecpos(:,3) = abs(elec.elecpos(:,3));
for i=1:nchan
  elec.elecpos(i,:) = 10*elec.elecpos(i,:)/norm(elec.elecpos(i,:));
end
elec.chanpos = elec.elecpos;
elec.tra = eye(nchan);

geom = [];
geom(1).pos = elec.elecpos*71/85; geom(1).unit = 'cm'; % brain
geom(2).pos = elec.elecpos*72/85; geom(2).unit = 'cm'; % csf
geom(3).pos = elec.elecpos*79/85; geom(3).unit = 'cm'; % skull
geom(4).pos = elec.elecpos*85/85; geom(4).unit = 'cm'; % scalp
headmodel = ft_headmodel_concentricspheres(geom, 'conductivity', [0.33 1.00 0.042 0.33]);
end


function headmodel = local_nestedspheres()
% four nested triangulated spheres as a BEM-style head model, ordered outside-in, in mm.
% no solver is attached; this is used for the geometry-only inside test and as input to
% ft_headmodel_openmeeg for the gated fits.
radii = [90 82 80 78]; % scalp, skull, csf, brain (mm)
[upos, tri] = mesh_sphere(642, 'icosahedron');
bnd = struct('pos', {}, 'tri', {});
for i=1:numel(radii)
  bnd(i).pos = upos*radii(i);
  bnd(i).tri = tri;
  bnd(i).unit = 'mm';
end
headmodel = [];
headmodel.bnd  = bnd;
headmodel.cond = [0.43 0.01 1.79 0.33];
headmodel.type = 'openmeeg';
headmodel.unit = 'mm';
end


function elec = local_projectelec(elec, radius)
% place electrodes on a sphere of the given radius (mm), keeping their directions, so they
% sit on the BEM scalp surface
dir = elec.elecpos ./ sqrt(sum(elec.elecpos.^2, 2));
elec.elecpos = dir*radius;
elec.chanpos = elec.elecpos;
elec.unit = 'mm';
end


function comp = local_topo2comp(topo, elec)
% wrap one or more spatial topographies as a FieldTrip component structure, so that
% ft_dipolefitting can fit one source per column with cfg.component
ncomp = size(topo,2);
comp = [];
comp.label    = arrayfun(@(k) sprintf('comp%03d', k), (1:ncomp)', 'UniformOutput', false);
comp.topolabel = elec.label;
comp.topo     = topo;       % nchan x ncomp
comp.unmixing = pinv(topo); % ncomp x nchan
comp.time     = {0};
comp.trial    = {zeros(ncomp,1)};
comp.elec     = elec;
end
