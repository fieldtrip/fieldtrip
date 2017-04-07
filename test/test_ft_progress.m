function test_ft_progress

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_progress

% testing parameters
pause_delay = 0.01;
num_iterations = 42;

ft_progress('init', 'etf',     'Please wait...');
for i=1:num_iterations
  ft_progress(i/num_iterations, 'Processing event %d from %d', i, num_iterations);
  pause(pause_delay);
end
ft_progress('close');

ft_progress('init', 'gui',     'Please wait...');
for i=1:num_iterations
  ft_progress(i/num_iterations, 'Processing event %d from %d', i, num_iterations);
  pause(pause_delay);
end
ft_progress('close');

ft_progress('init', 'text',     'Please wait...');
for i=1:num_iterations
  ft_progress(i/num_iterations, 'Processing event %d from %d', i, num_iterations);
  pause(pause_delay);
end
ft_progress('close');

ft_progress('init', 'dial',     'Please wait...');
for i=1:num_iterations
  ft_progress(i/num_iterations, 'Processing event %d from %d', i, num_iterations);
  pause(pause_delay);
end
ft_progress('close');

ft_progress('init', 'none',     'Please wait...');
for i=1:num_iterations
  ft_progress(i/num_iterations, 'Processing event %d from %d', i, num_iterations);
  pause(pause_delay);
end
ft_progress('close');

ft_progress('init', 'textbar',     'Please wait...');
for i=1:num_iterations
  ft_progress(i/num_iterations, 'Processing event %d from %d', i, num_iterations);
  pause(pause_delay);
end
ft_progress('close');

end
