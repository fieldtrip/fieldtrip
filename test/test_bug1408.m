function test_bug1408

% MEM 4gb
% WALLTIME 00:10:00

% TEST ft_preprocessing preproc ft_preproc_bandpassfilter ft_preproc_bandstopfilter ft_preproc_dftfilter ft_preproc_highpassfilter ft_preproc_lowpassfilter

nchans   = 200;
nsamples = 1e6;

dat = rand(nchans,nsamples); datmean = rand(nchans,1);
tic
dat = dat - datmean(:,ones(1,nsamples));
t(1) = toc;


dat = rand(nchans,nsamples); datmean = rand(nchans,1);
tic
dat = dat - repmat(datmean,[1 nsamples]);
t(2) = toc;

dat = rand(nchans,nsamples); datmean = rand(nchans,1);
tic
for ichan = 1:nchans
  dat(ichan,:) = dat(ichan,:) - datmean(ichan);
end
t(3) = toc;

dat = rand(nchans,nsamples); datmean = rand(nchans,1);
tic
dat = dat';
for ichan = 1:nchans
  dat(:,ichan) = dat(:,ichan) - datmean(ichan);
end
dat = dat';
t(4) = toc;

dat = rand(nchans,nsamples); datmean = rand(nchans,1);
tic
dat = bsxfun(@minus, dat, datmean);
t(5) = toc;

dat = rand(nchans,nsamples); datmean = rand(nchans,1);
tic
for isample = 1:nsamples
  dat(:,isample) = dat(:, isample) - datmean;
end
t(6) = toc;

[minval, minindx] = min(t);
if minindx~=6
  warning('unexpected winner of the speed test');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dat = rand(10,300)+ 30;
Fs  = 1000;
Fl  = 50;
Flp = 10;
Fhp = 70;
Fbp = [40 70];

filt = ft_preproc_bandpassfilter(dat,Fs,Fbp); assert(all(mean(filt,2)<1)); % the DC component should disappear
filt = ft_preproc_bandstopfilter(dat,Fs,Fbp); assert(all(mean(filt,2)>0));
filt = ft_preproc_dftfilter(dat,Fs,Fl);       assert(all(mean(filt,2)>0));
filt = ft_preproc_highpassfilter(dat,Fs,Fhp); assert(all(mean(filt,2)<1)); % the DC component should disappear
filt = ft_preproc_lowpassfilter(dat,Fs,Flp);  assert(all(mean(filt,2)>0));
