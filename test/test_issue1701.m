function test_issue1701

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY ft_debug ft_info ft_notice ft_warning ft_notification

global ft_default
ft_default.notification = [];

n = 5000;

%%

logging = @fprintf;

stopwatch = tic;
for i=1:n
  logging('hello 1\n');
end
t1 = toc(stopwatch);

%%

% logging = @ft_debug;
logging = @ft_info;
% logging = @ft_notice;
% logging = @ft_warning;

stopwatch = tic;
for i=1:n
  logging('hello 2\n');
end
t2 = toc(stopwatch);

%%

disp(t2/t1);

% fail if it gets more than than 200x slower
% this may be due to an inadverted slowdown of a lower-level function such as FT_VERSION 
assert(t2/t1<200);
