function test_issue1618

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_volumewrite

%%
mri      = [];
mri.dim  = [2 2 2];
mri.transform = eye(4);
mri.unit = 'mm';
mri.dat  = randn(mri.dim);

datatype = {'bit1', 'uint8', 'int16', 'int32', 'float', 'double'};
fname    = fullfile(tempdir, 'issue1618');

cfg           = [];
cfg.parameter = 'dat';

for k = 1:numel(datatype)
  cfg.datatype = datatype{k};
  cfg.filename = [fname, '_', datatype{k}];
  ft_volumewrite(cfg, mri);
end

