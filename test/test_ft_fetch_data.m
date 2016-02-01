function test_ft_fetch_data

% MEM 4gb
% WALLTIME 00:10:00

% TEST ft_fetch_data

% this function primarily tests the speed of ft_fetch_data

nchan     = 300;
fsample   = 1200;
ntrials   = 1000;

cachedata.label = {};
for c=1:nchan
  cachedata.label{c, 1}     = num2str(c);
end
cachedata.fsample        = fsample;
cachedata.time           = {};
cachedata.trial          = {};
cachedata.cfg            = [];
cachedata.sampleinfo = zeros(0,2);

endsample = 0;
for t=1:ntrials
  % variable trial length
  nsamples  = 2*ceil(1000000*rand/fsample);
  nsamples  = 1200;
  begsample = endsample + ceil(10*rand);
  endsample = begsample + nsamples-1;
  
  % make dummy data
  cachedata.sampleinfo(end+1,:) = [begsample endsample];
  cachedata.trial{end+1} = rand(nchan, nsamples);
  cachedata.time{end+1} = (1:nsamples)/cachedata.fsample;  
end

fprintf('avg. #samples: %.2f +/- %.2f\n', mean(cachedata.sampleinfo(:, 2) - cachedata.sampleinfo(:, 1)), std(cachedata.sampleinfo(:, 2) - cachedata.sampleinfo(:, 1)));
maxsamples = endsample + 1000;
% profile clear;
% profile on;
tic
for i = 1:100
  begsample = ceil(rand*maxsamples);
  endsample = begsample+ceil(4*rand*fsample);
  tmp = ft_fetch_data(cachedata, 'begsample', begsample', 'endsample', endsample);
end
toc
% profile report;


%% check with only one trial

nchan     = 300;
fsample   = 1200;
ntrials   = 1;

cachedata.label = {};
for c=1:nchan
  cachedata.label{c, 1}     = num2str(c);
end
cachedata.fsample        = fsample;
cachedata.time           = {};
cachedata.trial          = {};
cachedata.cfg            = [];
cachedata.sampleinfo = zeros(0,2);
nsamples  = 1200;

begsample = 100;
endsample = 100+1200-1;
  
% make dummy data
cachedata.sampleinfo(end+1,:) = [begsample endsample];
cachedata.trial{end+1} = rand(nchan, nsamples);
cachedata.time{end+1} = (1:nsamples)/cachedata.fsample;  

% completely outside, all should be nans
begsample = 40;
endsample = 80;
tmp = ft_fetch_data(cachedata, 'begsample', begsample', 'endsample', endsample);
if any(~isnan(tmp(:)))
  error('not all nans!'); 
end
% only end in
begsample = 90;
endsample = 100;
tmp = ft_fetch_data(cachedata, 'begsample', begsample', 'endsample', endsample);
if any(isnan(tmp(:, end))) || any(any(~isnan(tmp(:, 1:end-1))))
  error('nans are wrong!'); 
end

% all in
begsample = 100;
endsample = 101;
tmp = ft_fetch_data(cachedata, 'begsample', begsample', 'endsample', endsample);
if any(isnan(tmp(:)))
  error('nans are there!!'); 
end

% only beginning in
begsample = 100+1200-1;
endsample = 100+1200;
tmp = ft_fetch_data(cachedata, 'begsample', begsample', 'endsample', endsample);
if any(isnan(tmp(:, 1))) || any(any(~isnan(tmp(:, 2:end))))
  error('nans are wrong!'); 
end

% all out
begsample = 100+1200;
endsample = 100+1200+1;
tmp = ft_fetch_data(cachedata, 'begsample', begsample', 'endsample', endsample);
if any(~isnan(tmp(:)))
  error('not all nans!'); 
end


end
