function test_ft_progress

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_progress
% DATA no

% testing parameters
num_iterations = 100;
pause_delay = 5/num_iterations;

%%
fprintf('------------------ BEGIN none -------------------------------\n')
ft_progress('init', 'none',     'Please wait...');
for i=1:num_iterations
  ft_progress(i/num_iterations, 'Processing event %d from %d', i, num_iterations);
  pause(pause_delay);
end
ft_progress('close');
fprintf('------------------  END  -------------------------------\n')

%%
fprintf('------------------ BEGIN text -------------------------------\n')
ft_progress('init', 'text',     'Please wait...');
for i=1:num_iterations
  ft_progress(i/num_iterations, 'Processing event %d from %d', i, num_iterations);
  pause(pause_delay);
end
ft_progress('close');
fprintf('------------------  END  -------------------------------\n')

%%
fprintf('------------------ BEGIN textbar -------------------------------\n')
ft_progress('init', 'textbar',     'Please wait...');
for i=1:num_iterations
  ft_progress(i/num_iterations, 'Processing event %d from %d', i, num_iterations);
  pause(pause_delay);
end
ft_progress('close');
fprintf('------------------  END  -------------------------------\n')

%%
fprintf('------------------ BEGIN etf -------------------------------\n')
ft_progress('init', 'etf',     'Please wait...');
for i=1:num_iterations
  ft_progress(i/num_iterations, 'Processing event %d from %d', i, num_iterations);
  pause(pause_delay);
end
ft_progress('close');
fprintf('------------------  END  -------------------------------\n')

%%
fprintf('------------------ BEGIN dial -------------------------------\n')
ft_progress('init', 'dial',     'Please wait...');
for i=1:num_iterations
  ft_progress(i/num_iterations, 'Processing event %d from %d', i, num_iterations);
  pause(pause_delay);
end
ft_progress('close');
fprintf('------------------  END  -------------------------------\n')

%%
fprintf('------------------ BEGIN gui -------------------------------\n')
ft_progress('init', 'gui',     'Please wait...');
for i=1:num_iterations
  ft_progress(i/num_iterations, 'Processing event %d from %d', i, num_iterations);
  pause(pause_delay);
end
ft_progress('close');
fprintf('------------------  END  -------------------------------\n')

