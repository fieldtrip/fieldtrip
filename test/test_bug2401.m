function test_bug2401

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_prepare_headmodel ft_prepare_vol_sens ft_read_vol read_ctf_hdm

hdr = ft_read_header(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds'));
sens = hdr.grad;

chansel1 = 1:184;
chansel2 = 1:151;
chansel3 = 151:-1:1;
chansel4 = [27 130 100 33 25 89 147 20 44 26 28 23 63 58 101 38 75 41 51 140 46 21];

cfg = [];
cfg.method = 'localspheres';
cfg.hdmfile = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds/default.hdm');
cfg.channel = hdr.label(chansel1);
vol = ft_prepare_headmodel(cfg); % this basically reads it from file

% an alternative approach would be to use
% hdm = ft_read_vol(cfg.hdmfile);

[vol1p, sens1p] = ft_prepare_vol_sens(vol, sens, 'channel', sens.label(chansel1));
[vol2p, sens2p] = ft_prepare_vol_sens(vol, sens, 'channel', sens.label(chansel2));
[vol3p, sens3p] = ft_prepare_vol_sens(vol, sens, 'channel', sens.label(chansel3));
[vol4p, sens4p] = ft_prepare_vol_sens(vol, sens, 'channel', sens.label(chansel4));

% channel MLF41 is number 25 in the original grad structure
sel1 = find(strcmp('MLF41', sens1p.label));
sel2 = find(strcmp('MLF41', sens2p.label));
sel3 = find(strcmp('MLF41', sens3p.label));
sel4 = find(strcmp('MLF41', sens4p.label));

r1 = vol1p.r(sens1p.tra(sel1,:)~=0,:);
r2 = vol2p.r(sens2p.tra(sel2,:)~=0,:);
r3 = vol3p.r(sens3p.tra(sel3,:)~=0,:);
r4 = vol4p.r(sens4p.tra(sel4,:)~=0,:);

o1 = vol1p.o(sens1p.tra(sel1,:)~=0,:);
o2 = vol2p.o(sens2p.tra(sel2,:)~=0,:);
o3 = vol3p.o(sens3p.tra(sel3,:)~=0,:);
o4 = vol4p.o(sens4p.tra(sel4,:)~=0,:);

assert(isequal(r1, r2));
assert(isequal(r1, r3));
assert(isequal(r1, r4));

assert(isequal(o1, o2));
assert(isequal(o1, o3));
assert(isequal(o1, o4));
