function test_pull401

% WALLTIME 00:10:00
% MEM 1gb

level = {
  'emergency'
  'alert'
  'critical'
  'error'
  'warning'
  'warning_once'
  'notice'
  'info'
  'debug'};

cfg = [];


fprintf('======================== using ft_notice ======================== \n');
for i=1:numel(level)
  fprintf('LEVEL = %s\n', level{i});
  cfg.feedback = level{i};
  ft_notice('this is a notice\n');
end

fprintf('======================== using ft_info ======================== \n');
for i=1:numel(level)
  fprintf('LEVEL = %s\n', level{i});
  cfg.feedback = level{i};
  ft_info('this is an info message\n');
end


%%
fprintf('======================== using ft_warning ======================== \n');

% please note when running/edditing this test script that it refers to fieldtrip/test/private/ft_warning
% and not to fieldtrip/utilities/ft_warning
ft_warning('-clear')

stopwatch = tic;
while toc(stopwatch)<5
  disp(toc(stopwatch))
  ft_warning('this is a warning', 1);
  drawnow
  pause(0.1);
end
