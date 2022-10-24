function test_ft_denoise_hfc(testdata)

% WALLTIME XX:XX:XX
% MEM 8gb
% DEPENDENCY test_ft_denoise_hfc

switch nargin
    case 0
        locate_function = which('test_ft_denoise_hfc');
        % load data
        load(fullfile(locate_function,'..','..','..','project_zeus',...
            'data_UCL_OPM.mat'));
    case 1
        % If user inputs a filename, use path specific for DCCN
        % Find this data: https://doi.org/10.17605/OSF.IO/CJNXH
        load(dccnpath('...data_UCL_OPM.mat'));
end


% Test of Harmonic Field Correction (HFC) where L = 1
cfg               = [];
cfg.L             = 1;
cfg.residualcheck = 'yes';
[data_out]        = ft_denoise_hfc(cfg,data);

% Check L = 2
cfg               = [];
cfg.L             = 2;
cfg.residualcheck = 'yes';
[data_out]        = ft_denoise_hfc(cfg,data);

% Check L = 3
cfg               = [];
cfg.L             = 3;
cfg.residualcheck = 'yes';
[data_out]        = ft_denoise_hfc(cfg,data);

% Check L = 1 without residual check
cfg               = [];
cfg.L             = 1;
cfg.residualcheck = 'no';
[data_out]        = ft_denoise_hfc(cfg,data);




