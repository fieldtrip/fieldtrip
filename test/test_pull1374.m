function test_pull1374

% WALLTIME 00:15:00
% MEM 4gb
% DEPENDENCY xsens_mvnx motion_c3d
% DATA private

%%

filename=dccnpath('/project/3031000.02/test/original/motion/xsens/pull1374.mvnx');

ft_read_header(filename);
ft_read_data(filename);

cfg = [];
cfg.dataset = filename;
data = ft_preprocessing(cfg);

%%

if ft_platform_supports('matlabversion', '2018b', inf)
  
  filename=dccnpath('/project/3031000.02/test/original/motion/xsens/sub-03_rec-01_take-02_motion.c3d');
  
  ft_read_header(filename);
  ft_read_data(filename);
  
  cfg = [];
  cfg.dataset = filename;
  data = ft_preprocessing(cfg);
  
else
  
  warning('the ezc3d mex file is not compatible with older MATLAB versions');
  
end