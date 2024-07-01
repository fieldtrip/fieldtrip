function test_issue2082

% MEM 1gb
% WALLTIME 00:30:00
% DEPENDENCY ft_checkconfig
% DATA public

global ft_default

cfg = [];
cfg.dataset = string(dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds'));
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');

filename = string(dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.shape'));
hs = ft_read_headshape(filename);

% create a dummy cfg
cfg.a = 'a';
cfg.b = "b";
cfg.c = 1;
cfg.d.e = "e";
cfg.d.f = 1;
cfg.d.g.h = 'h';
cfg.d.g.i = "i";
cfg.d.g.j = {"a" "b"};
k.a = "a";
k.b = 'b';
cfg.d.g.k = {k 1}; % this for now does not work

cfg.checkstring = 'no';
cfgout1 = ft_checkconfig(cfg);
c = struct2cell(cfgout1);
s = fieldnames(cfg);
t = cellfun(@class, c, 'UniformOutput', false);
assert(any(strcmp(t, 'string')));

cfg.checkstring = 'yes';
cfgout2 = ft_checkconfig(cfg);

c = struct2cell(cfgout2);
t = cellfun(@class, c, 'UniformOutput', false);
assert(~any(strcmp(t, 'string')));

c = struct2cell(cfgout2.d);
t = cellfun(@class, c, 'UniformOutput', false);
assert(~any(strcmp(t, 'string')));
