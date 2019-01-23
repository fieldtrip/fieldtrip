function test_pull832

% WALLTIME 00:10:00
% MEM 3gb

% TEST ft_read_sens ft_read_header mne2grad

%%

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/pull832'));

%%

filename = 'sub-01_ses-meg_task-facerecognition_run-01_meg.fif';

hdrx = ft_read_header(filename, 'coilaccuracy', []);
hdr0 = ft_read_header(filename, 'coilaccuracy', 0);
hdr1 = ft_read_header(filename, 'coilaccuracy', 1);
hdr2 = ft_read_header(filename, 'coilaccuracy', 2);

%%

gradx = ft_read_sens(filename, 'senstype', 'meg', 'coilaccuracy', []);
grad0 = ft_read_sens(filename, 'senstype', 'meg', 'coilaccuracy', 0);
grad1 = ft_read_sens(filename, 'senstype', 'meg', 'coilaccuracy', 1);
grad2 = ft_read_sens(filename, 'senstype', 'meg', 'coilaccuracy', 2);

if false
  % this only ran once on 12 Oct 2018 with the old code
  gradx_old = ft_read_sens(filename, 'senstype', 'meg', 'coilaccuracy', []);
  grad0_old = ft_read_sens(filename, 'senstype', 'meg', 'coilaccuracy', 0);
  grad1_old = ft_read_sens(filename, 'senstype', 'meg', 'coilaccuracy', 1);
  grad2_old = ft_read_sens(filename, 'senstype', 'meg', 'coilaccuracy', 2);
  
  save gradx_old gradx_old
  save grad0_old grad0_old
  save grad1_old grad1_old
  save grad2_old grad2_old
else
  load gradx_old
  load grad0_old
  load grad1_old
  load grad2_old
end

assert( isequal(gradx, gradx_old));
assert(~isequal(grad0, grad0_old));
assert(~isequal(grad1, grad1_old));
assert(~isequal(grad2, grad2_old));

%%

vol = [];
vol.o = [0 0 4];
vol.r = 10;
vol.unit = 'cm';

dippos = [0 0 9];
dipmom = [0 1 0]';

lfx = ft_compute_leadfield(dippos, gradx, vol) * dipmom;

%%

vol = [];
vol.o = [0 0 0.04];
vol.r = 0.10;
vol.unit = 'm';

dippos = [0 0 0.09];
dipmom = [0 1 0]';

lf0 = ft_compute_leadfield(dippos, grad0, vol) * dipmom;
lf1 = ft_compute_leadfield(dippos, grad1, vol) * dipmom;
lf2 = ft_compute_leadfield(dippos, grad2, vol) * dipmom;

%%

lf0_old = ft_compute_leadfield(dippos, grad0_old, vol) * dipmom;
lf1_old = ft_compute_leadfield(dippos, grad1_old, vol) * dipmom;
lf2_old = ft_compute_leadfield(dippos, grad2_old, vol) * dipmom;

%%

figure
ft_plot_dipole(dippos, dipmom, 'unit', 'm')
ft_plot_sens(grad0)

%%

figure
plot([lf0, lf0_old]);
hold on
plot(lf0-lf0_old, 'k');

relvar = (norm(lf0-lf0_old) ./ norm(lf0))^2
assert(relvar<0.03);

figure
plot([lf1, lf1_old]);
hold on
plot(lf1-lf1_old, 'k');

relvar = (norm(lf1-lf1_old) ./ norm(lf1))^2
assert(relvar<0.03);

figure
plot([lf2, lf2_old]);
hold on
plot(lf2-lf2_old, 'k');

relvar = (norm(lf2-lf2_old) ./ norm(lf2))^2
assert(relvar<0.03);

%%

pla1 = (0:101)*3 + 1;
pla2 = (0:101)*3 + 2;
pla = [pla1 pla2];
mag = (0:101)*3 + 3;

relvar = (norm(lf0(mag)-lf0_old(mag)) ./ norm(lf0(mag)))^2
relvar = (norm(lf0(pla)-lf0_old(pla)) ./ norm(lf0(pla)))^2

