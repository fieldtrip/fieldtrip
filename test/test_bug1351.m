function test_bug1351

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug1315
% TEST ft_topoplotER ft_topoplotTFR ft_checksize ft_checkconfig

% load the test data of just two subjects
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1351.mat'));

% plot topographically - note: plot will look strange
cfg = [];
cfg.layout = 'CTF275.lay';

max_it = 10;
figure;

profile on -history -timer real;
t = inf(max_it, 1);
for i=1:max_it
tic;
evalc('ft_topoplotER(cfg, tmpdata)'); % surpress output
t(i) = toc;
end
profile off;

stats = profile('info');
time = [stats.FunctionTable.TotalTime];
[tmp, idx] = sort(time);

fprintf('\n');
fprintf('-------------------------------------------------------------\n');
fprintf('-------------- Functions and computation time ---------------\n');
fprintf('-------------------------------------------------------------\n');
fprintf('\n');
fprintf('Function\t\t\tTime\t\t\tRelativeTime\n');
for i=0:15
  fprintf('%s\t\t\t%.2fs\t\t\t%.2f%%\n', ...
    stats.FunctionTable(idx(end-i)).FunctionName, ...
    stats.FunctionTable(idx(end-i)).TotalTime, ...
    100*...
      stats.FunctionTable(idx(end-i)).TotalTime./...
      stats.FunctionTable(idx(end)).TotalTime);
end
fprintf('\n');


% fprintf('Plotting took on average %.2fms (%.2fms to %.2fms)...\n', 1000*mean(t), 1000*min(t), 1000*max(t));
% close all;
% if min(t) > 5 % note that usually we plot >2 subjects and >15 channels
%   error('This is unacceptable')
% else
%   fprintf('This is fine!\n');
% end

