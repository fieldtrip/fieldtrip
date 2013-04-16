<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> enhancement - started making a test script for bug 1967
=======
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967
function test_bug1967
=======
% function test_bug1967
>>>>>>> enhancement - started making a test script for bug 1967
<<<<<<< HEAD
<<<<<<< HEAD
=======
function test_bug1967
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967
=======
>>>>>>> enhancement - started making a test script for bug 1967
=======
=======
function test_bug1967
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967
=======
function test_bug1967
>>>>>>> merge conflict
=======
% function test_bug1967
>>>>>>> enhancement - started making a test script for bug 1967
=======
function test_bug1967
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967

% TEST test_bug1967
% TEST ft_prepare_vol_sens

%% construct a segmentation
nx = 101;
ny = 101;
nz = 101;

seg = zeros(nx, ny, nz);
[posx, posy, posz] = ndgrid(1:nx, 1:ny,1:nz);
posx = posx-(nx-1)/2-1;
posy = posy-(nx-1)/2-1;
posz = posz-(nx-1)/2-1;

distance = sqrt(posx.^2 + posy.^2 + posz.^2);
seg(distance<45) = 1;
seg(distance<40) = 2;
seg(distance<30) = 3;

volume.dim = [nx ny nz];
volume.seg = seg;
volume.seglabel = {'skin' 'skull' 'brain'};
volume.transform = [
  1 0 0 -51 
  0 1 0 -51
  0 0 1 -51 
  0 0 0 1
  ];

cfg = [];
cfg.funparameter = 'seg';
ft_sourceplot(cfg, volume);

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967
%% convert it into a mesh 
=======
%% convert it into a head model
>>>>>>> enhancement - started making a test script for bug 1967
=======
%% convert it into a mesh 
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967
<<<<<<< HEAD
=======
%% convert it into a mesh 
=======
%% convert it into a head model
>>>>>>> enhancement - started making a test script for bug 1967
>>>>>>> enhancement - started making a test script for bug 1967
=======
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967
=======
%% convert it into a mesh 
>>>>>>> merge conflict
=======
%% convert it into a head model
>>>>>>> enhancement - started making a test script for bug 1967
=======
%% convert it into a mesh 
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967
cfg = [];
cfg.tissue = {'skin' 'skull' 'brain'};
cfg.method = 'hexahedral';
mesh = ft_prepare_mesh(cfg, volume);

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967
%% convert it into a head model
cfg = [];
cfg.tissue = {'skin' 'skull' 'brain'};
cfg.conductivity = [1 1/80 1];
cfg.method = 'simbio';
headmodel = ft_prepare_headmodel(cfg, mesh);

%% make some electrodes
[pnt, tri] = icosahedron42;
pnt = pnt(pnt(:,3)>0,:);
pnt = pnt*55; % not precisely fitting on the mesh
elec = [];
elec.chanpos = pnt;
elec.elecpos = pnt;
for i=1:size(pnt,1)
  elec.label{i} = num2str(i);
end
elec.unit = 'mm';

figure
ft_plot_mesh(mesh, 'edgeonly', 1)
ft_plot_sens(elec)

%% prepare the volume conductor and electrodes for leadfield computation
[vol, sens] = ft_prepare_vol_sens(headmodel, elec)

% elec is the original one, sens is the one after projecting
if isequal(sens.elecpos, elec.elecpos)
  error('the electrodes were not projected');
end
=======
=======
<<<<<<< HEAD
=======
>>>>>>> enhancement - started making a test script for bug 1967
=======
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967
=======
>>>>>>> merge conflict
%% convert it into a head model
cfg = [];
cfg.tissue = {'skin' 'skull' 'brain'};
cfg.conductivity = [1 1/80 1];
cfg.method = 'simbio';
headmodel = ft_prepare_headmodel(cfg, mesh);

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> enhancement - started making a test script for bug 1967
=======
%% make some electrodes
[pnt, tri] = icosahedron42;
pnt = pnt(pnt(:,3)>0,:);
pnt = pnt*55; % not precisely fitting on the mesh
elec = [];
elec.chanpos = pnt;
elec.elecpos = pnt;
for i=1:size(pnt,1)
  elec.label{i} = num2str(i);
end
elec.unit = 'mm';

figure
ft_plot_mesh(mesh, 'edgeonly', 1)
ft_plot_sens(elec)

%% prepare the volume conductor and electrodes for leadfield computation
[vol, sens] = ft_prepare_vol_sens(headmodel, elec)

% elec is the original one, sens is the one after projecting
if isequal(sens.elecpos, elec.elecpos)
  error('the electrodes were not projected');
end
<<<<<<< HEAD
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967
=======
=======
=======
>>>>>>> enhancement - started making a test script for bug 1967
=======
%% convert it into a head model
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967
cfg = [];
cfg.tissue = {'skin' 'skull' 'brain'};
cfg.conductivity = [1 1/80 1];
cfg.method = 'simbio';
headmodel = ft_prepare_headmodel(cfg, mesh);

<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> enhancement - started making a test script for bug 1967
=======
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967
>>>>>>> enhancement - started making a test script for bug 1967
=======
=======
>>>>>>> merge conflict
=======
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967
%% make some electrodes
[pnt, tri] = icosahedron42;
pnt = pnt(pnt(:,3)>0,:);
pnt = pnt*55; % not precisely fitting on the mesh
elec = [];
elec.chanpos = pnt;
elec.elecpos = pnt;
for i=1:size(pnt,1)
  elec.label{i} = num2str(i);
end
elec.unit = 'mm';

figure
ft_plot_mesh(mesh, 'edgeonly', 1)
ft_plot_sens(elec)

%% prepare the volume conductor and electrodes for leadfield computation
[vol, sens] = ft_prepare_vol_sens(headmodel, elec)

% elec is the original one, sens is the one after projecting
if isequal(sens.elecpos, elec.elecpos)
  error('the electrodes were not projected');
end
<<<<<<< HEAD
=======
>>>>>>> enhancement - started making a test script for bug 1967
=======
>>>>>>> enhancement - extended test script for http://bugzilla.fcdonders.nl/show_bug.cgi?id=1967

