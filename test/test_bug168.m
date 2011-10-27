function test_bug168

% TEST test_bug168
% TEST ft_realtime_topography

if isempty(which('ft_realtime_topography'))
  addpath(fullfile(fileparts(which('ft_defaults')), 'realtime/example'));
end

cfg = [];
cfg.dataset = '/home/common/matlab/fieldtrip/data/Subject01.ds';
cfg.bufferdata = 'first';
cfg.layout = 'CTF151.lay';

ft_realtime_topography(cfg);


