function inspect_ft_sourceplot_interactive

% MEM 8gb
% WALLTIME 00:10:00
% DEPENDENCY ft_sourceplot_interactive ft_plot_mesh_interactive

%% load first dataset: results from MNE tutorial

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/minimumnormestimate/source.mat'));

%%

cfg = [];
ft_sourceplot_interactive(cfg, sourceFC, sourceFIC);

fprintf('Inspect the generated plot (all defaults); press a key to continue...\n');
pause;
close all;


%%

cfg = [];
cfg.parameter = 'pow';
cfg.operation = 'subtract';
difference = ft_math(cfg, sourceFC, sourceFIC);

cfg = [];
cfg.data_labels = {'Congruent', 'Incongruent', 'Difference'};
cfg.has_diff = true;
ft_sourceplot_interactive(cfg, sourceFC, sourceFIC, difference);

fprintf('Inspect the generated plot (labels present; has difference); press a key to continue...\n');
pause;
close all;

%% warp data to MNI space (to test atlas)

% transform is approximate, but good enough for test purposes
transform_mni2ctf = [0.0173    0.9945    0.1038   36.6736
                    -0.9995    0.0144    0.0284    0.0389
                     0.0268   -0.1042    0.9942   52.4980
                          0         0         0    1.0000];
       
transf = pinv(transform_mni2ctf);
                         
sourceFC_mni = ft_transform_geometry(transf, sourceFC);
sourceFIC_mni = ft_transform_geometry(transf, sourceFIC);
difference_mni = ft_transform_geometry(transf, difference);

cfg = [];
cfg.data_labels = {'Congruent', 'Incongruent', 'Difference'};
cfg.has_diff = true;
[ftver, ftpath] = ft_version();
cfg.atlas = fullfile(ftpath, 'template', 'atlas', 'aal', 'ROI_MNI_V4.nii');
ft_sourceplot_interactive(cfg, sourceFC_mni, sourceFIC_mni, difference_mni);

fprintf('Inspect the generated plot (additionally also atlas label present); press a key to continue...\n');
pause;
close all;

end
