function test_bug1288

% MEM 1500mb
% WALLTIME 00:10:00

% this function serves to create planar gradient data and combined planar
% gradient data, for testing purposes of a fix for bug 1288

load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat'));

cfg = [];
cfg.method = 'distance';
cfg.grad   = data.grad;
neighbours = ft_prepare_neighbours(cfg);

cfg = [];
cfg.planarmethod = 'sincos';
cfg.neighbours   = neighbours;
datap = ft_megplanar(cfg, data);

cfg = [];
datac = ft_combineplanar(cfg, datap);

cfg = [];
cfg.grad = datac.grad;
lay = ft_prepare_layout(cfg);

% this appears to run through without error, but
cfg = [];
ft_layoutplot(cfg, datac);
% results in a figure that has channels ML* (left) towards the nose and MR* (right)
% towards the back of the head. So 90 degrees rotated clockwise.

if false
  % this is the correct reference solution
  cfg = [];
  cfg.layout = 'CTF151.lay';
  lay = ft_prepare_layout(cfg);
end

posL = lay.pos(strcmp(lay.label, 'MLT43'),:);
posR = lay.pos(strcmp(lay.label, 'MRT43'),:);

dx = posR(1)-posL(1);
dy = posR(2)-posL(2);

assert(dx>0.5); % should be large
assert(dy<0.1); % should be small

% the problem arose because ft_prepare_layout does not detect that it is a CTF system and consequently it does not set the default rotation to 90 degrees
% so the actual underlying checks that need to be performed are
assert(ft_senstype(data.grad, 'ctf151'))
assert(ft_senstype(datap.grad, 'ctf151_planar'))
assert(ft_senstype(datac.grad, 'ctf151')); % this assertion failed in revision 5725 and probably all previous versions




