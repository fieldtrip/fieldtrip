function test_ft_getopt

% TEST: ft_getopt

success = true;

success = success && isequal(ft_getopt({'key1', 'val1'}, 'key1'), 'val1')
success = success && isequal(ft_getopt({'key1', 'val1'}, 'key2'), [])
success = success && isequal(ft_getopt({'key1', 'val1'}, 'key2', 'default'), 'default')
success = success && isequal(ft_getopt({'key1', 'val1'}, 'key2', 'default'), 'default')
clear cfg
cfg.key1 = 'val1';
success = success && isequal(ft_getopt(cfg, 'key1'), 'val1')
success = success && isequal(ft_getopt(cfg, 'key2'), [])
success = success && isequal(ft_getopt(cfg, 'key2', 'default'), 'default')
success = success && isequal(ft_getopt(cfg, 'key2', 'default'), 'default')

success = success && isequal(ft_getopt({'key1', []}, 'key1'), [])
success = success && isequal(ft_getopt({'key1', []}, 'key1', 'default'), 'default')
clear cfg
cfg.key1 = [];
success = success && isequal(ft_getopt(cfg, 'key1'), [])
success = success && isequal(ft_getopt(cfg, 'key1', 'default'), 'default')

success = success && isequal(ft_getopt({'key1', 'val1', 'key2', 'val2', 'key3', 'val3'}, 'key3'), 'val3')
success = success && isequal(ft_getopt({'key1', 'val1', 'key2', 'val2', 'key3', 'val3'}, 'key4'), [])
success = success && isequal(ft_getopt({'key1', 'val1', 'key2', 'val2', 'key3', 'val3'}, 'key4', 'default'), 'default')
clear cfg
cfg.key1 = 'val1';
cfg.key2 = 'val2';
cfg.key3 = 'val3';
success = success && isequal(ft_getopt(cfg, 'key3'), 'val3')
success = success && isequal(ft_getopt(cfg, 'key4'), [])
success = success && isequal(ft_getopt(cfg, 'key4', 'default'), 'default')

if ~success
  error('there is a problem with ft_getopt');
end

