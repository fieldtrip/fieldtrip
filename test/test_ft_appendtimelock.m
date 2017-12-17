function test_ft_appendtimelock

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_appendtimelock

% Johanna Zumer

% make some dummy timelock structures
tlock.label = {'1';'2'};
tlock.time  = 1:5;
tlock.dimord = 'rpt_chan_time';
tlock.trial = randn(10,2,5);

% concat over trials
tlock2=tlock;
tlock2.trial=randn(8,2,5);

cfg = [];
tlockapp = ft_appendtimelock(cfg,tlock,tlock2);
assert(all( size(tlockapp.trial)==[18 2 5]))

cfg = [];
cfg.channel={'1'};
tlockapp = ft_appendtimelock(cfg,tlock,tlock2);
assert(all( size(tlockapp.trial)==[18 1 5]))

% concat over channels
clear tlock2
tlock1.label=tlock.label(1);
tlock2.label=tlock.label(2);
tlock1.trial=tlock.trial(:,1,:);
tlock2.trial=tlock.trial(:,2,:);
tlock1.time=tlock.time;
tlock2.time=tlock.time;
tlock1.dimord=tlock.dimord;
tlock2.dimord=tlock.dimord;

cfg = [];
tlockapp = ft_appendtimelock(cfg,tlock,tlock2);
assert(all( size(tlockapp.trial)==[20 1 5]))

cfg = [];
tlockapp = ft_appendtimelock(cfg,tlock1,tlock2);
assert(all( size(tlockapp.trial)==[10 2 5]))

%%  Test using other combinations of tlock from real data
clear tlock*
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/eventrelatedaveraging/dataFC_LP.mat'));

cfg=[];
tlock1=ft_timelockanalysis(cfg,dataFC_LP);

cfg=[];
cfg.channel={'MLT*'};
tlockl=ft_timelockanalysis(cfg,dataFC_LP);
cfg=[];
cfg.channel={'MRT*'};
tlockr=ft_timelockanalysis(cfg,dataFC_LP);

% should concat over channels
cfg=[];
tlockapp=ft_appendtimelock(cfg,tlockl,tlockr);

try
  cfg=[];
  cfg.appenddim='rpt';
  tlockapp=ft_appendtimelock(cfg,tlockl,tlockr);
  catchflag=0;
catch
  catchflag=1;
end
assert(logical(catchflag))

% should use only subset of 'l' channels and concat over avg
cfg=[];
tlockapp=ft_appendtimelock(cfg,tlockl,tlock1);

try
  cfg=[];
  cfg.appenddim='chan';
  tlockapp=ft_appendtimelock(cfg,tlockl,tlock1);
  catchflag=0;
catch
  catchflag=1;
end
assert(logical(catchflag))

% new tlock with keeptrials
cfg=[];
cfg.keeptrials = 'yes';
tlock1=ft_timelockanalysis(cfg,dataFC_LP);
cfg.channel={'MLT*'};
tlockl=ft_timelockanalysis(cfg,dataFC_LP);
cfg.channel={'MRT*'};
tlockr=ft_timelockanalysis(cfg,dataFC_LP);

% should concat over channels
cfg=[];
tlockapp=ft_appendtimelock(cfg,tlockl,tlockr);

try
  cfg=[];
  cfg.appenddim='rpt';
  tlockapp=ft_appendtimelock(cfg,tlockl,tlockr);
  catchflag=0;
catch
  catchflag=1;
end
assert(logical(catchflag))

% should use only subset of 'l' channels and concat over avg
cfg=[];
tlockapp=ft_appendtimelock(cfg,tlockl,tlock1);


% now test for numerical inaccurracies, should concatenate across 'rpt'
tlock2          = tlock1;
tlock2.time     = tlock1.time+0.0000001;
cfg=[];
tlockapp       = ft_appendtimelock(cfg, tlock1, tlock2);

%% test for data with labels shuffled around

tlock2=tlock1;
tlock2.label = tlock1.label(randperm(length(tlock1.label)));

tlock2.trial = randn(size(tlock1.trial));

cfg=[];
tlockshuffled = ft_appendtimelock(cfg, tlock1, tlock2);


% now check whether channels were correctly appended
[a,b] = match_str(tlock1.label, tlockshuffled.label);
x1 = tlock1.trial(:,a(1),:);
y1 = tlockshuffled.trial(1:size(x1,1),b(1),:);

[a,b] = match_str(tlock2.label, tlockshuffled.label);
x2 = tlock2.trial(:,a(1),:);
y2 = tlockshuffled.trial([size(x1,1)+1]:end,b(1),:);

if ~all(x1(:) == y1(:)) || ~all(x2(:) == y2(:))
  error('data was wrongly appended when channel labels are differently ordered in input arguments');
end

%% testing with covariance subfield
% FIXME: create test with .cov 
