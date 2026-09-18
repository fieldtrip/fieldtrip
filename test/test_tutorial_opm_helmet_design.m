function test_tutorial_opm_helmet_design(datadir)

% MEM 2gb
% WALLTIME 00:20:00
% DEPENDENCY ft_read_headshape ft_electroderealign ft_defacemesh
% DATA public

% see http://www.fieldtriptoolbox.org/tutorial/preproc/denoising_opm
% this test function corresponds to the version on the wiki at 10 September
% 2026, NOTE that it has been stripped down to the parts that can run
% non-interactively, so that it is a part of the daily test batch. Also the
% function can be used for quick testing (e.g. just prior to a toolkit) to
% check whether the basic functionality works, provided the version of the
% tutorial online has not deviated much since the current m-file has been
% created.

if nargin<1
  datadir = dccnpath('/project/3031000.02/external/download/tutorial/opm_helmet_design');
end
pwdir = pwd;
cd(datadir);

%
headshape = ft_read_headshape('spherical-head.stl');
helmet = ft_read_headshape('spherical-helmet.stl');

figure
ft_plot_mesh(headshape, 'facecolor', 'skin_light', 'axes', true)
ft_plot_mesh(helmet, 'facecolor', [0.7 0.7 0.7]) % grey
ft_headlight % or lighting phong; camlight

%headshape.coordsys = 'ctf'; 
helmet.coordsys = 'ctf'; 

%
nas = [+100 0 0];
ini = [-100 0 0];
lpa = [0 +100 0];
rpa = [0 -100 0];

%
headshape.fid.pos = [
    nas
    ini
    lpa
    rpa
];

headshape.fid.label = {
    'nas'
    'ini'
    'lpa'
    'rpa'
};

figure
ft_plot_headshape(headshape, 'facecolor', 'skin_light', 'axes', true)
alpha 0.5 % slightly transparent
ft_headlight

%
cfg = [];
cfg.fiducial.nas = nas;
cfg.fiducial.ini = ini;
cfg.fiducial.lpa = lpa;
cfg.fiducial.rpa = rpa;
cfg.method = '1020';
cfg.feedback = 'yes';
headshape.coordsys = 'als'; % added by JM to avoid interactive stuff
elec = ft_electrodeplacement(cfg, headshape);

%
chansel = ft_channelselection({'eeg1020', '-Fpz', '-Oz'}, elec.label);

cfg = [];
cfg.elec = elec;
cfg.channel = chansel;
cfg.outwardshift = 10/2; % the helmet is 10 mm thick, the bottom of the sensor holder will be halfway in

cfg.template = 'fieldline_sensor.stl';
[tmpcfg, sensor] = ft_sensorplacement(cfg, headshape);

cfg.template = 'fieldline_holder.stl';
[tmpcfg, holder] = ft_sensorplacement(cfg, headshape);

cfg.template = 'fieldline_hole.stl';
[tmpcfg, hole] = ft_sensorplacement(cfg, headshape);

cfg.template = 'fieldline_padding.stl';
[tmpcfg, padding] = ft_sensorplacement(cfg, headshape);

%
figure
ft_plot_headshape(headshape, 'facecolor', 'skin', 'facealpha', 1, 'axes', true);
ft_plot_mesh(helmet, 'facecolor', 'lightgray', 'facealpha', 0.5, 'edgecolor', 'none');
ft_plot_mesh(sensor, 'facecolor', 'r', 'facealpha', 1, 'edgecolor', 'none');
ft_plot_mesh(holder, 'facecolor', 'g', 'facealpha', 1, 'edgecolor', 'none');
% ft_plot_mesh(hole, 'facecolor', 'b', 'facealpha', 0.5, 'edgecolor', 'none');
% ft_plot_mesh(padding, 'facecolor', 'm', 'facealpha', 0.5, 'edgecolor', 'none');
ft_headlight

pwdir = pwd;
cd(tempdir);
mkdir('spherical');
for i=1:numel(sensor)
    disp(chansel{i});
    filename = sprintf('spherical/sensor-%s.stl', chansel{i});
    ft_write_headshape(filename, sensor(i), 'fileformat', 'stl');
    filename = sprintf('spherical/holder-%s.stl', chansel{i});
    ft_write_headshape(filename, holder(i), 'fileformat', 'stl');
    filename = sprintf('spherical/hole-%s.stl', chansel{i});
    ft_write_headshape(filename, hole(i), 'fileformat', 'stl');
    filename = sprintf('spherical/padding-%s.stl', chansel{i});
    ft_write_headshape(filename, padding(i), 'fileformat', 'stl');
end
cd(pwdir);

%
headshape = ft_read_headshape('flattenedspherical-head.stl');
helmet = ft_read_headshape('flattenedspherical-helmet.stl');

%
nas = [+100 0 -10];
ini = [-100 0 -10];
lpa = [ 0 +80 -10];
rpa = [ 0 -80 -10];

headshape.fid.pos = [
    nas
    ini
    lpa
    rpa
];

headshape.fid.label = {
    'nas'
    'ini'
    'lpa'
    'rpa'
};

figure
ft_plot_headshape(headshape, 'facecolor', 'skin_light', 'axes', true)
alpha 0.5 % slightly transparent
ft_headlight

%
headshape.coordsys = 'als';
helmet.coordsys = 'als';

%
cfg = [];
cfg.fiducial.nas = nas;
cfg.fiducial.ini = ini;
cfg.fiducial.lpa = lpa;
cfg.fiducial.rpa = rpa;
cfg.method = '1020';
cfg.feedback = 'yes';
elec = ft_electrodeplacement(cfg, headshape);

%
chansel = ft_channelselection({'eeg1020', '-Fpz', '-Oz'}, elec.label);

cfg = [];
cfg.elec = elec;
cfg.channel = chansel;
cfg.outwardshift = 10/2; % the helmet is 10 mm thick, the bottom of the sensor holder will be halfway in

cfg.template = 'fieldline_sensor.stl';
[tmpcfg, sensor] = ft_sensorplacement(cfg, headshape);

cfg.template = 'fieldline_holder.stl';
[tmpcfg, holder] = ft_sensorplacement(cfg, headshape);

cfg.template = 'fieldline_hole.stl';
[tmpcfg, hole] = ft_sensorplacement(cfg, headshape);

cfg.template = 'fieldline_padding.stl';
[tmpcfg, hole] = ft_sensorplacement(cfg, headshape);

%
figure
ft_plot_headshape(headshape, 'facecolor', 'skin', 'facealpha', 0.5, 'axes', 1);
ft_plot_mesh(helmet, 'facecolor', 'lightgray', 'facealpha', 0.5, 'edgecolor', 'none');
ft_plot_mesh(sensor, 'facecolor', 'r', 'facealpha', 1, 'edgecolor', 'none');
ft_plot_mesh(holder, 'facecolor', 'g', 'facealpha', 1, 'edgecolor', 'none');
% ft_plot_mesh(hole, 'facecolor', 'b', 'facealpha', 0.5, 'edgecolor', 'none');
% ft_plot_mesh(padding, 'facecolor', 'm', 'facealpha', 0.5, 'edgecolor', 'none');
ft_headlight

% % for i=1:numel(sensor)
% %     disp(chansel{i});
% %     filename = sprintf('flattenedspherical-sensor-%s.stl', chansel{i});
% %     ft_write_headshape(filename, sensor(i), 'fileformat', 'stl');
% %     filename = sprintf('flattenedspherical-holder-%s.stl', chansel{i});
% %     ft_write_headshape(filename, holder(i), 'fileformat', 'stl');
% %     filename = sprintf('flattenedspherical-hole-%s.stl', chansel{i});
% %     ft_write_headshape(filename, hole(i), 'fileformat', 'stl');
% %     filename = sprintf('flattenedspherical-padding-%s.stl', chansel{i});
% %     ft_write_headshape(filename, padding(i), 'fileformat', 'stl');
% % end

%
mri = ft_read_mri('individual.nii');

ft_determine_coordsys(mri, 'interactive', 'no')
rotate3d

%
cfg = [];
cfg.method = 'ortho';
cfg.flip = 'no'; % important for identifying voxel indices
ft_sourceplot(cfg, mri);

nas_vox = [ 102 217 123 ];
ini_vox = [ 90   18 104 ];
lpa_vox = [ 20  135 103 ];
rpa_vox = [ 173 121  95 ];

%
cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = 'ctf';
cfg.fiducial.nas = nas_vox;
cfg.fiducial.lpa = lpa_vox;
cfg.fiducial.rpa = rpa_vox;
mri_realigned = ft_volumerealign(cfg, mri);

%
nas = ft_transform_geometry(mri.transform, nas_vox);
ini = ft_transform_geometry(mri.transform, ini_vox);
lpa = ft_transform_geometry(mri.transform, lpa_vox);
rpa = ft_transform_geometry(mri.transform, rpa_vox);

cfg = [];
cfg.method = 'ortho';
cfg.location = nas; % nas, ini, lpa, rpa
cfg.locationcoordinates = 'head';
ft_sourceplot(cfg, mri_realigned);

%
cfg = [];
mri_resliced = ft_volumereslice(cfg, mri_realigned);

cfg = [];
cfg.method = 'ortho';
cfg.location = [0 0 0];
cfg.locationcoordinates = 'head';
ft_sourceplot(cfg, mri_resliced);

%
disp(mri_resliced.cfg)

cfg = [];
cfg.xrange = [ -97.5000 157.5000] - 30;
cfg.yrange = [-127.5000 127.5000];
cfg.zrange = [ -87.5000 167.5000] - 5;
mri_resliced = ft_volumereslice(cfg, mri_realigned);

cfg = [];
cfg.method = 'ortho';
cfg.location = [0 0 0];
cfg.locationcoordinates = 'head';
ft_sourceplot(cfg, mri_resliced);

%
cfg = [];
cfg.output = 'scalp';
mri_segmented = ft_volumesegment(cfg, mri_resliced);

cfg = [];
cfg.method = 'projectmesh';
cfg.numvertices = 4000;
headshape = ft_prepare_mesh(cfg, mri_segmented);

%
headshape.fid.pos = [
    nas
    lpa
    rpa
    ini
];

headshape.fid.label = {
    'nas'
    'lpa'
    'rpa'
    'ini'
};

figure
ft_plot_headshape(headshape, 'facecolor', 'skin', 'facealpha', 0.7, 'fidmarker', 'o', 'fidcolor', 'k', 'fidlabel', true, 'fidsize', 24)
ft_headlight

%
cfg = [];
cfg.funparameter = 'scalp';
ft_sourceplot(cfg, mri_segmented)

%
mri_segmented.scalp(:,:,1:50) = 0;
mri_segmented.airgap = imdilate(mri_segmented.scalp,  strel('sphere', 1));
mri_segmented.helmet = imdilate(mri_segmented.airgap, strel('sphere', 5));

%
mri_indexed = ft_checkdata(mri_segmented, 'segmentationstyle', 'indexed');

cfg = [];
cfg.funparameter = 'tissue';
cfg.location = [0 0 0];
cfg.locationcoordinates = 'head';
cfg.atlas = mri_indexed;
ft_sourceplot(cfg, mri_indexed)

%
cfg = [];
cfg.method = 'projectmesh';
cfg.numvertices = 4000;

cfg.tissue = 'airgap'; % this makes a surface from the outside of the "airgap" part
tmp = removefields(mri_segmented, {'scalp', 'helmet'}); % FIXME this is a hack that should be resolved
inside = ft_prepare_mesh(cfg, tmp);

cfg.tissue = 'helmet';
tmp = removefields(mri_segmented, {'scalp', 'airgap'}); % FIXME this is a hack that should be resolved
outside = ft_prepare_mesh(cfg, tmp);

%
figure
ft_plot_headshape(headshape, 'facecolor', 'skin', 'facealpha', 0.5, 'edgecolor', 'none');
ft_plot_mesh(inside, 'facecolor', 'r', 'facealpha', 0.5, 'edgecolor', 'none');
ft_headlight
rotate3d

%
figure
ft_plot_headshape(headshape, 'facecolor', 'skin', 'facealpha', 0.5, 'edgecolor', 'none');
ft_plot_mesh(outside, 'facecolor', 'g', 'facealpha', 0.5, 'edgecolor', 'none');
ft_headlight
rotate3d

%
cfg = [];
cfg.fiducial.nas = nas;
cfg.fiducial.ini = ini;
cfg.fiducial.lpa = lpa;
cfg.fiducial.rpa = rpa;
cfg.method = '1020';
elec = ft_electrodeplacement(cfg, headshape);

%
chansel = ft_channelselection({'eeg1020', '-Fpz', '-Oz'}, elec.label);

cfg = [];
cfg.elec = elec;
cfg.channel = chansel;
cfg.outwardshift = 5/2 + 1; % the helmet is 5 mm thick, plus the 1 mm airgap

cfg.template = 'fieldline_sensor.stl';
[tmpcfg, sensor] = ft_sensorplacement(cfg, headshape);

cfg.template = 'fieldline_holder.stl';
[tmpcfg, holder] = ft_sensorplacement(cfg, headshape);

cfg.template = 'fieldline_hole.stl';
[tmpcfg, hole] = ft_sensorplacement(cfg, headshape);

cfg.template = 'fieldline_padding.stl';
[tmpcfg, hole] = ft_sensorplacement(cfg, headshape);

%
figure
ft_plot_headshape(headshape, 'facecolor', 'skin', 'facealpha', 1, 'axes', 1);
ft_plot_mesh(outside, 'facecolor', 'lightgray', 'facealpha', 0.5, 'edgecolor', 'none');
ft_plot_mesh(sensor, 'facecolor', 'r', 'facealpha', 1, 'edgecolor', 'none');
ft_plot_mesh(holder, 'facecolor', 'g', 'facealpha', 1, 'edgecolor', 'none');
% ft_plot_mesh(hole, 'facecolor', 'b', 'facealpha', 0.5, 'edgecolor', 'none');
% ft_plot_mesh(padding, 'facecolor', 'm', 'facealpha', 0.5, 'edgecolor', 'none');
ft_headlight

%
pwdir = pwd;
cd(tempdir);
mkdir('individual');
ft_write_headshape('individual/helmet-inside.stl', inside, 'fileformat', 'stl');
ft_write_headshape('individual/helmet-outside.stl', outside, 'fileformat', 'stl');

for i=1:numel(sensor)
disp(chansel{i});
    filename = sprintf('individual/sensor-%s.stl', chansel{i});
    ft_write_headshape(filename, sensor(i), 'fileformat', 'stl');
    filename = sprintf('individual/holder-%s.stl', chansel{i});
    ft_write_headshape(filename, holder(i), 'fileformat', 'stl');
    filename = sprintf('individual/hole-%s.stl', chansel{i});
    ft_write_headshape(filename, hole(i), 'fileformat', 'stl');
    filename = sprintf('individual/padding-%s.stl', chansel{i});
    ft_write_headshape(filename, padding(i), 'fileformat', 'stl');
end
cd(pwdir);

%
nsubj = 10;

mri = cell(1,nsubj);
fiducial = cell(1,nsubj);

for i=1:nsubj
    fprintf('------------------------------- %d -------------------------------\n', i);

    mrifile = sprintf('population/subject%03d.nii', i);
    mri{i} = ft_read_mri(mrifile);
    mri{i}.coordsys = 'ctf';
    
    fidfile = sprintf('population/fiducial%03d.mat', i);
    fiducial{i} = load(fidfile);
end

%
ft_determine_coordsys(mri{1}, 'interactive', 'no')

%
cfg = [];
cfg.method = 'ortho';
cfg.location = fiducial{3}.nas;
cfg.locationcoordinates = 'head';
ft_sourceplot(cfg, mri{1});

%
mri_resliced = {};
mri_segmented = {};

for i=1:nsubj
    fprintf('------------------------------- %d -------------------------------\n', i);
    cfg = [];
    cfg.xrange = [ -97.5000 157.5000] - 30;
    cfg.yrange = [-127.5000 127.5000];
    cfg.zrange = [ -87.5000 167.5000] - 5;
    mri_resliced{i} = ft_volumereslice(cfg, mri{i});

    cfg = [];
    cfg.output = 'scalp';
    mri_segmented{i} = ft_volumesegment(cfg, mri_resliced{i});
end

%
cfg = [];
cfg.method = 'ortho';
cfg.location = fiducial{1}.nas;
cfg.funparameter = 'scalp';
ft_sourceplot(cfg, mri_segmented{1});

%
mri_averaged = rmfield(mri_segmented{1}, 'cfg');
for i=2:nsubj
    mri_averaged.scalp = mri_averaged.scalp + mri_segmented{i}.scalp;
end
mri_averaged.scalp = mri_averaged.scalp/nsubj;

%
cfg = [];
cfg.method = 'ortho';
cfg.location = [0 0 0];
cfg.funparameter = 'scalp';
ft_sourceplot(cfg, mri_averaged);

%
mri_90percentile       = mri_averaged;
mri_90percentile.scalp = mri_averaged.scalp>=0.1;

%
cfg = [];
cfg.method = 'ortho';
cfg.location = [0 0 0];
cfg.funparameter = 'scalp';
ft_sourceplot(cfg, mri_90percentile);

% make it perfectly symmetric
mri_90percentile.scalp =  mri_90percentile.scalp | flip(mri_90percentile.scalp, 2);

% we now do the same as for the individual MRI
mri_segmented = mri_90percentile;

% we still need to make the headshape mesh, so don't modify the scalp segmentation directly
tmp = mri_segmented.scalp;
tmp(:,:,1:50) = 0;

mri_segmented.airgap = imdilate(tmp,                  strel('sphere', 1));
mri_segmented.helmet = imdilate(mri_segmented.airgap, strel('sphere', 5));

%
cfg = [];
cfg.method = 'projectmesh';
cfg.numvertices = 4000;

cfg.tissue = 'scalp';
headshape = ft_prepare_mesh(cfg, mri_segmented);

cfg.tissue = 'airgap'; % this makes a surface from the outside of the "airgap" part
tmp = removefields(mri_segmented, {'scalp', 'helmet'}); % FIXME this is a hack that should be resolved
inside = ft_prepare_mesh(cfg, tmp);

cfg.tissue = 'helmet';
tmp = removefields(mri_segmented, {'scalp', 'airgap'}); % FIXME this is a hack that should be resolved
outside = ft_prepare_mesh(cfg, tmp);

%
figure
ft_plot_axes([], 'unit', 'mm', 'coordsys', 'als');
for i=1:nsubj
    ft_plot_mesh(fiducial{i}.nas, 'vertexcolor', 'k');
    ft_plot_mesh(fiducial{i}.lpa, 'vertexcolor', 'k');
    ft_plot_mesh(fiducial{i}.rpa, 'vertexcolor', 'k');
    ft_plot_mesh(fiducial{i}.ini, 'vertexcolor', 'k');
end

%
nas_avg = fiducial{1}.nas;
ini_avg = fiducial{1}.ini;
lpa_avg = fiducial{1}.lpa;
rpa_avg = fiducial{1}.rpa;

for i=2:nsubj
    nas_avg = nas_avg + fiducial{i}.nas;
    ini_avg = ini_avg + fiducial{i}.ini;
    lpa_avg = lpa_avg + fiducial{i}.lpa;
    rpa_avg = rpa_avg + fiducial{i}.rpa;
end

nas_avg = nas_avg/nsubj;
ini_avg = ini_avg/nsubj;
lpa_avg = lpa_avg/nsubj;
rpa_avg = rpa_avg/nsubj;

% average the position of left and right ear
tmp1 = (lpa_avg + [1 -1 1] .* rpa_avg)/2;
tmp2 = (rpa_avg + [1 -1 1] .* lpa_avg)/2;
% now we can safely overwrite them
lpa_avg = tmp1;
rpa_avg = tmp2;

% the nose and inion should be exactly on the y=0 plane
nas_avg(2) = 0;
ini_avg(2) = 0;

%
headshape.fid.pos = [
    nas_avg
    ini_avg
    lpa_avg
    rpa_avg
    ];

headshape.fid.label = {
    'nas'
    'ini'
    'lpa'
    'rpa'
    };

figure
ft_plot_headshape(headshape, 'axes', 'on', 'facecolor', 'skin', 'facealpha', 0.5, 'fidmarker', '.', 'fidcolor', 'k', 'fidlabel', true, 'fidsize', 24)
ft_headlight

%
cfg = [];
cfg.fiducial.nas = nas_avg + [+10 0 0];
cfg.fiducial.ini = ini_avg + [-10 0 0];
cfg.fiducial.lpa = lpa_avg + [0 +10 0];
cfg.fiducial.rpa = rpa_avg + [0 -10 0];
cfg.method = '1020';
elec = ft_electrodeplacement(cfg, headshape);

%
chansel = ft_channelselection({'eeg1020', '-Fpz', '-Oz'}, elec.label);

%
chansel = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'T7', 'C3', 'Cz', 'C4', 'T8', 'P7', 'P3', 'Pz', 'P4', 'P8', 'O1', 'O2'};

%
cfg = [];
cfg.layout = 'EEG1005.lay';
cfg.channel = {'all'}; % modify this to your selection
cfg.skipcomnt = 'yes';
cfg.skipscale = 'yes';
cfg.feedback = 'yes';
layout = ft_prepare_layout(cfg);

%
headshape = ft_read_headshape('spherical-head.stl');
headshape.coordsys = 'ctf';

%
[headshape.tri, headshape.pos] = reducepatch(headshape.tri, headshape.pos, 800);

%
nas = [+100 0 0];
ini = [-100 0 0];
lpa = [0 +100 0];
rpa = [0 -100 0];

cfg = [];
cfg.fiducial.nas = nas;
cfg.fiducial.ini = ini;
cfg.fiducial.lpa = lpa;
cfg.fiducial.rpa = rpa;
cfg.method = 'equidistant';
cfg.numelec = 64;
cfg.maxiter = 500;
cfg.feedback = 'yes';
elec = ft_electrodeplacement(cfg, headshape);

%
sensor = ft_read_headshape('fieldline_sensor.stl');
holder = ft_read_headshape('fieldline_holder.stl');

figure
ft_plot_mesh(sensor, 'facecolor', 'r', 'edgecolor', 'none')
ft_plot_mesh(holder, 'facecolor', 'lightgray', 'edgecolor', 'none')
ft_plot_axes([], 'unit', 'mm')
ft_headlight

axis on; grid on
xlabel('x')
ylabel('y')
zlabel('z')

%
headshape = ft_read_headshape('spherical-head.stl');

nas = [+100 0 0];
ini = [-100 0 0];
lpa = [0 +100 0];
rpa = [0 -100 0];

cfg = [];
cfg.fiducial.nas = nas;
cfg.fiducial.ini = ini;
cfg.fiducial.lpa = lpa;
cfg.fiducial.rpa = rpa;
cfg.method = '1020';
cfg.feedback = 'no';
headshape.coordsys = 'als';
elec = ft_electrodeplacement(cfg, headshape);

chansel = ft_channelselection({'eeg1020', '-Fpz', '-Oz'}, elec.label);

cfg = [];
cfg.elec = elec;
cfg.channel = chansel;
cfg.outwardshift = 2; % two mm away from the surface 
cfg.template = 'fieldline_sensor.stl';
[outcfg, sensor] = ft_sensorplacement(cfg, headshape);

figure
ft_plot_headshape(headshape, 'facecolor', 'skin', 'facealpha', 1, 'axes', 1);
ft_plot_mesh(sensor, 'facecolor', 'r', 'facealpha', 1, 'edgecolor', 'none');
ft_headlight

%
cfg = [];
cfg.rotx = zeros(19,1);
cfg.roty = zeros(19,1);
cfg.rotz = zeros(19,1);
cfg.elec = elec;
cfg.channel = chansel;
cfg.outwardshift = 2; % two mm away from the surface 
cfg.template = 'fieldline_sensor.stl';
[tmpcfg, sensor] = ft_sensorplacement(cfg, headshape);

figure
ft_plot_headshape(headshape, 'facecolor', 'skin', 'facealpha', 1, 'axes', 1);
ft_plot_mesh(sensor, 'facecolor', 'r', 'facealpha', 1, 'edgecolor', 'none');
ft_headlight

%
cfg = [];
cfg.rotx = zeros(19,1);
cfg.roty = ones(19,1) * 45; % degrees
cfg.rotz = zeros(19,1);
cfg.elec = elec;
cfg.channel = chansel;
cfg.outwardshift = 2; % two mm away from the surface 
cfg.template = 'fieldline_sensor.stl';
[tmpcfg, sensor] = ft_sensorplacement(cfg, headshape);

figure
ft_plot_headshape(headshape, 'facecolor', 'skin', 'facealpha', 1, 'axes', 1);
ft_plot_mesh(sensor, 'facecolor', 'r', 'facealpha', 1, 'edgecolor', 'none');
ft_headlight

%
selT7 = find(strcmp(tmpcfg.channel, 'T7'));
selT8 = find(strcmp(tmpcfg.channel, 'T8'));

cfg = [];
cfg.rotx = outcfg.rotx; % copy this from the previous call 
cfg.roty = outcfg.roty;
cfg.rotz = outcfg.rotz;

% do not rotate these around the y-axis, which connects the ears
cfg.rotz(selT7) = 0;
cfg.rotz(selT8) = 0;
% rotate plusminus 90 degrees around the x-axis, which goes to the nose
cfg.rotx(selT7) = +90;
cfg.rotx(selT8) = -90;
% ... the remainder would be the same as before

%
for i=1:numel(elec.label)
    elec.elecori(i,:) = elec.elecpos(i,:) - [0 0 40];
    elec.elecori(i,:) = elec.elecori(i,:) / norm(elec.elecori(i,:)); % unit length
end

%
cfg = [];
cfg.rotx = zeros(19,1);
cfg.roty = zeros(19,1);
cfg.rotz = ones(19,1) * 45; % degrees
cfg.elec = elec;
cfg.channel = chansel;
cfg.outwardshift = 2; % two mm away from the surface 
cfg.template = 'fieldline_sensor.stl';
[tmpcfg, sensor] = ft_sensorplacement(cfg, headshape);

figure
ft_plot_headshape(headshape, 'facecolor', 'skin', 'facealpha', 1, 'axes', 1);
ft_plot_mesh(sensor, 'facecolor', 'r', 'facealpha', 1, 'edgecolor', 'none');
ft_headlight

%
correction = {
    'Cz'    0 0   0
    'T7'    0 0 -30
    'T8'    0 0 +30
    % ...
    };

cfg = [];
cfg.elec = elec;
cfg.channel = chansel;  % subset of 19 locations
cfg.rotx = outcfg.rotx; % 19x1 vector
cfg.roty = outcfg.roty; % 19x1 vector
cfg.rotz = outcfg.rotz; % 19x1 vector

% chansel is a cell-array with the selected locations of the 1020 system at which sensors are placed
for i=1:size(correction,1)
  lab = correction{i,1};
  dx  = correction{i,2};
  dy  = correction{i,3};
  dz  = correction{i,4};
  cfg.rotx(strcmp(chansel, lab)) = cfg.rotx(strcmp(chansel, lab)) + dx;
  cfg.roty(strcmp(chansel, lab)) = cfg.roty(strcmp(chansel, lab)) + dy;
  cfg.rotz(strcmp(chansel, lab)) = cfg.rotz(strcmp(chansel, lab)) + dz;
end

%
headshape = ft_read_headshape('spherical-head.stl');
headshape.coordsys = 'ctf';

nas = [+100 0 0];
ini = [-100 0 0];
lpa = [0 +100 0];
rpa = [0 -100 0];

cfg = [];
cfg.fiducial.nas = nas;
cfg.fiducial.ini = ini;
cfg.fiducial.lpa = lpa;
cfg.fiducial.rpa = rpa;
cfg.method = '1020';
cfg.feedback = 'no';
elec = ft_electrodeplacement(cfg, headshape);

chansel = ft_channelselection({'eeg1020', '-Fpz', '-Oz'}, elec.label);

cfg = [];
cfg.elec = elec;
cfg.channel = chansel;
cfg.template = 'fieldline_sensor.stl';
cfg.outwardshift = 2; % two mm away from the surface
[outcfg, sensor] = ft_sensorplacement(cfg, headshape);

%
opm_single = [];
opm_single.label = {
    'x'
    'y'
    'z'
};
opm_single.coilpos = [
    0 0 0
    0 0 0
    0 0 0
];
opm_single.coilori = [
    1 0 0
    0 1 0
    0 0 1
];
opm_single.tra = eye(3);

%
cfg = [];
cfg.elec = elec;
cfg.channel = chansel;
cfg.template = opm_single;
cfg.outwardshift = 2 + 1 + 5; % IMPORTANT see below
[outcfg, all_opm] = ft_sensorplacement(cfg, headshape);

%
lab = {};
grad = [];
grad.label = {};
grad.coilpos = zeros(0,3);
grad.coilori = zeros(0,3);
for i=1:length(chansel)
    lab{1} = [outcfg.channel{i} '_' all_opm(i).label{1}]; % _x
    lab{2} = [outcfg.channel{i} '_' all_opm(i).label{2}]; % _y
    lab{3} = [outcfg.channel{i} '_' all_opm(i).label{3}]; % _z
    grad.label = cat(1, grad.label, lab(:));
    grad.coilpos = cat(1, grad.coilpos, all_opm(i).coilpos);
    grad.coilori = cat(1, grad.coilori, all_opm(i).coilori);
end
grad.tra = eye(length(grad.label));

% We can plot the resulting OPM sensor positions for all channels with ft_plot_sens.
%
figure
ft_plot_headshape(headshape, 'facecolor', 'skin', 'facealpha', 1, 'axes', 1);
ft_plot_sens(grad);
ft_headlight

%
% or only with the z-oriented channels, including their labels.
%
ft_plot_sens(grad, 'chanindx', endsWith(grad.label, 'z'), 'label', 'label')

%
cfg = [];
cfg.grad = grad; % this contains x, y, and z channels
cfg.channel = grad.label(endsWith(grad.label, 'z'));
cfg.rotate = 90;
cfg.feedback = 'yes';
layout = ft_prepare_layout(cfg);

%
cfg = [];
cfg.grad = grad; % this contains x, y, and z channels
cfg.rotate = 90;
cfg.feedback = 'no';

cfg.channel = grad.label(endsWith(grad.label, 'x'));
layout_x = ft_prepare_layout(cfg);

cfg.channel = grad.label(endsWith(grad.label, 'y'));
layout_y = ft_prepare_layout(cfg);

cfg.channel = grad.label(endsWith(grad.label, 'z'));
layout_z = ft_prepare_layout(cfg);

cfg = [];
layout_xyz = ft_appendlayout(cfg, layout_x, layout_y, layout_z);

figure
ft_plot_layout(layout_xyz)


cd(pwdir);
