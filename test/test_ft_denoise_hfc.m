function test_ft_denoise_hfc(testfile)

% WALLTIME 00:10:00
% MEM 6gb
% DEPENDENCY test_ft_denoise_hfc
% DATA private

% switch nargin
%   case 0
%     locate_function = which('test_ft_denoise_hfc');
%     % load data
%     load(fullfile(locate_function,'..','..','..','project_zeus',...
%       'data_UCL_OPM.mat'));
%   case 1
%     % If user inputs a filename, use path specific for DCCN
%     % Find this data: https://doi.org/10.17605/OSF.IO/CJNXH
%     load(dccnpath('/home/common/matlab/fieldtrip/data/test/pull2123/data_UCL_OPM.mat'));
% end
if nargin == 0
  testfile = dccnpath('/home/common/matlab/fieldtrip/data/test/pull2123/data_UCL_OPM.mat');
else
  % testfile is required in the input, as a full path to a mat-file
end
load(testfile);

% Test of Harmonic Field Correction (HFC) where L = 1
cfg               = [];
cfg.order         = 1;
cfg.residualcheck = 'yes';
[data_out]        = ft_denoise_hfc(cfg,data);

% Check L = 2
cfg               = [];
cfg.order         = 2;
cfg.residualcheck = 'yes';
[data_out]        = ft_denoise_hfc(cfg,data);

% Check L = 3
cfg               = [];
cfg.order         = 3;
cfg.residualcheck = 'yes';
[data_out]        = ft_denoise_hfc(cfg,data);

% Check L = 1 without residual check
cfg               = [];
cfg.order         = 1;
cfg.residualcheck = 'no';
[data_out]        = ft_denoise_hfc(cfg,data);




