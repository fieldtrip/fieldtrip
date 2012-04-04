function test_bug1408

% TEST test_bug1408

nchans   = 200;
nsamples = 1e6;

dat = rand(nchans,nsamples); datmean = rand(nchans,1);
tic
dat = dat - datmean(:,ones(1,nsamples));
t(1) = toc
 

dat = rand(nchans,nsamples); datmean = rand(nchans,1);
tic
dat = dat - repmat(datmean,[1 nsamples]);
t(2) = toc

dat = rand(nchans,nsamples); datmean = rand(nchans,1);
tic
for ichan = 1:nchans
  dat(ichan,:) = dat(ichan,:) - datmean(ichan);
end
t(3) = toc

dat = rand(nchans,nsamples); datmean = rand(nchans,1);
tic
dat = dat';
for ichan = 1:nchans
  dat(:,ichan) = dat(:,ichan) - datmean(ichan);
end
dat = dat';
t(4) = toc

dat = rand(nchans,nsamples); datmean = rand(nchans,1);
tic
dat = bsxfun(@minus, dat, datmean);
t(5) = toc

dat = rand(nchans,nsamples); datmean = rand(nchans,1);
tic
for isample = 1:nsamples
  dat(:,isample) = dat(:, isample) - datmean;
end
t(6) = toc

[minval, minindx] = min(t);
if minindx~=6
  error('unexpected winner of the speed test');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data          = [];
data.label    = {'1', '2', '3'};
data.trial{1} = rand(3,1000);
data.time{1}  = (0:999)/1000;
data.fsample  = 1000;

cfg = [];
cfg.demean = 'yes';
dataout = ft_preprocessing(cfg, data);

cfg = [];
cfg.detrend = 'yes';
dataout = ft_preprocessing(cfg, data);

cfg = [];
cfg.polyremoval = 'yes';
cfg.polyorder = 3;
dataout = ft_preprocessing(cfg, data);

