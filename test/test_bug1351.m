function test_bug1351
% TEST test_bug1315
% TEST ft_topoplotER ft_topoplotTFR ft_checksize checksize checkconfig
% TEST ft_checkconfig


% load test data of just two subjects
load test_bug1351

% plot topographically - note: plot will look strange
cfg = [];
cfg.layout = 'CTF275.lay';

max_it = 10;
figure;

t = inf(max_it, 1);
for i=1:max_it
tic;
evalc('ft_topoplotER(cfg, tmpdata)'); % surpress output
t(i) = toc;
end

fprintf('Plotting took on average %.2fms (%.2fms to %.2fms)...\n', 1000*mean(t), 1000*min(t), 1000*max(t));
close all;
if min(t) > 5 % note that usually we plot >2 subjects and >15 channels
  error('This is unacceptable')
else
  fprintf('This is fine!\n');
end