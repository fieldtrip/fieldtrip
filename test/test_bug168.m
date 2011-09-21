function test_bug168

% TEST test_bug168
% TEST ft_realtime_topography

cfg = [];
cfg.dataset = '/Users/robert/Manzana/data/MEG/Subject01.ds';
cfg.bufferdata = 'first';
cfg.layout = 'CTF151.lay';

ft_realtime_topography(cfg);


