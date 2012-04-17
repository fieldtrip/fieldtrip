function test_ft_progress

% TEST test_ft_progress
% TEST ft_progress

ft_progress('init', 'etf',     'Please wait...');
for i=1:42
  ft_progress(i/42, 'Processing event %d from %d', i, 42);
  pause(0.01);
end
ft_progress('close');

ft_progress('init', 'gui',     'Please wait...');
for i=1:42
  ft_progress(i/42, 'Processing event %d from %d', i, 42);
  pause(0.01);
end
ft_progress('close');

ft_progress('init', 'text',     'Please wait...');
for i=1:42
  ft_progress(i/42, 'Processing event %d from %d', i, 42);
  pause(0.01);
end
ft_progress('close');

ft_progress('init', 'dial',     'Please wait...');
for i=1:42
  ft_progress(i/42, 'Processing event %d from %d', i, 42);
  pause(0.01);
end
ft_progress('close');

ft_progress('init', 'none',     'Please wait...');
for i=1:42
  ft_progress(i/42, 'Processing event %d from %d', i, 42);
  pause(0.01);
end
ft_progress('close');

ft_progress('init', 'textbar',     'Please wait...');
for i=1:42
  ft_progress(i/42, 'Processing event %d from %d', i, 42);
  pause(0.01);
end
ft_progress('close');

end