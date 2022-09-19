function test_bug1800

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_defaults ft_selectdata ft_topoplotER

% this was reported by Giorgos Michalareas
%
% The data spans from roughly -4 to 4 sec.
% When the entire time of the data is used for the plot then the ft_topoplotER
% does not seem to perform any matching between the data channels and the layout
% channels and it seems to assume that the are already matched. This results in a
% wrong topology plot. 
%
% http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1800

% Load the timelocked data
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1800.mat'));

cfg = [];
cfg.avgovertime = 'yes';
tmpavg2=ft_selectdata(cfg,tmpavg1);  % average timelocked data across all time points

% Error  case - Plot the topoplot for timelocked data
cfg=[];
ft_topoplotER(cfg,tmpavg1);

% Error  case - Plot the topoplot for timelocked data averaged all time points
cfg=[];
ft_topoplotER(cfg,tmpavg2);

% Error  case - Plot the topoplot for timelocked data by setting xlim equal to the time range of the data;
cfg=[];
cfg.xlim=[tmpavg1.time(1) tmpavg1.time(end)];
ft_topoplotER(cfg,tmpavg1);

% Correct  case - Plot the topoplot for timelocked data by setting xlim to a SUBINTERVAL of the time range of the data;
cfg=[];
cfg.xlim=[tmpavg1.time(2) tmpavg1.time(end-1)];
ft_topoplotER(cfg,tmpavg1);

