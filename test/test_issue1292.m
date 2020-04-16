function test_issue1292

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_fetch_data

%% This one has identical overlap

data1 = [];
data1.label = {'1'};
data1.trial{1} = 1*ones(1,100);
data1.trial{2} = 1*ones(1,100);
data1.time{1} = 1:100;
data1.time{2} = 1:100;
data1.sampleinfo = [
  1 100
  91 190
  ];

%% This one has non-identical overlap

data2 = [];
data2.label = {'1'};
data2.trial{1} = 1*ones(1,100);
data2.trial{2} = 2*ones(1,100); % NOTE HERE
data2.time{1} = 1:100;
data2.time{2} = 1:100;
data2.sampleinfo = [
  1 100
  91 190
  ];

%%

dat1a = ft_fetch_data(data1, 'begsample', 1, 'endsample', 190);

%%

dat1b = ft_fetch_data(data1, 'begsample', 1, 'endsample', 190, 'allowoverlap', true);

%%

try
  haserror = true; % this will be set when it fails
  dat2a    = ft_fetch_data(data2, 'begsample', 1, 'endsample', 190);
  haserror = false; % this will be set when it passed
end

if ~haserror
  error('this should give an error');
end

%%

dat2b = ft_fetch_data(data2, 'begsample', 1, 'endsample', 190, 'allowoverlap', true);
