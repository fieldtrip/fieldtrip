function test_bug840

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug800
% TEST ft_preproc_resample

nchans = 13;
nsamples = 10000;

dat = randn(nchans, nsamples);

typ = {
  'double'
  'single'
  'int64'
  'int32'
  'int16'
  'int8'
  };

for i=1:length(typ)
  disp(typ{i})
  dat0 = cast(dat, typ{i});
  dat1 = ft_preproc_resample(dat0, 1000, 500, 'resample');
  dat2 = ft_preproc_resample(dat0, 1000, 500, 'decimate');
  dat3 = ft_preproc_resample(dat0, 1000, 500, 'downsample');
  
  % first column should be nchans
  assert(size(dat1,1)==nchans);
  assert(size(dat2,1)==nchans);
  assert(size(dat3,1)==nchans);
  
  % output type should remain constant
  assert(strcmp(class(dat1), typ{i}));
  assert(strcmp(class(dat2), typ{i}));
  assert(strcmp(class(dat3), typ{i}));
end


