function test_bug2800

% TEST project_elec

% WALLTIME 00:10:00
% MEM 500mb

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2800.mat');
load(filename);

% we need to cd into the private directory
privatedir = dccnpath('/home/common/matlab/fieldtrip/private');
cd(privatedir);

[el,prj] = project_elec(elec.elecpos, vol.bnd(1).pnt, vol.bnd(1).tri);
d        = sqrt(sum((elec.elecpos-prj).^2,2));

% take as cutoff 3 cm, i.e. 30 mm, this is a bit of a heuristic but is
% based on JM's inspection of the to-be-expected values in the correct case
assert(all(d<30));

