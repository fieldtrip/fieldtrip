function test_bug1529

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_plot_ortho

[ftver, ftpath] = ft_version;

%% Load the SPMtemplate 
mritem = ft_read_mri([ftpath '/external/spm8/templates/T1.nii']);
% mrichecked = ft_determine_coordsys(mritem) ;
mritem.coordsys = 'spm'; % define the coordinate system (why is this not part of ft_read_mri?)

%% Display the Anatomy via plot_ortho (defaults)
figure('name','MRI TEM');
ft_plot_ortho(mritem.anatomy);

%% Display the Anatomy via plot_ortho (at a certain location )
location_voxel = [46 64 37];
location_str   = sprintf('x: %.0f y: %.0f z:%.0f',location_voxel(1:3));
figure('name',['MRI TEM VOX' location_str ]);
ft_plot_ortho(mritem.anatomy,'location',location_voxel(1:3));

%% Display the Anatomy via plot_ortho (at the same certain location, using external transform )
location_mri = [0 0 0];
location_voxel = mritem.transform\[location_mri 1]';
location_str   = sprintf('x: %.0f y: %.0f z:%.0f',location_voxel(1:3));
figure('name',['MRI TEM MRI EXTERN ' location_str]);
ft_plot_ortho(mritem.anatomy,'location',location_voxel(1:3));

%% Display the Anatomy via plot_ortho (at the same certain location, using internal transform )
location_mri = [0 0 0];
location_str   = sprintf('x: %.0f y: %.0f z:%.0f',location_mri(1:3));
figure('name',['MRI TEM MRI INTERN ' location_str]);
ft_plot_ortho(mritem.anatomy,'location',location_mri(1:3),'transform',mritem.transform);

%% Display the Anatomy via plot_ortho (at the same certain location, using inverse internal transform )
location_mri = [0 0 0];
location_str   = sprintf('x: %.0f y: %.0f z:%.0f',location_mri(1:3));
figure('name',['MRI TEM MRI INTERN ' location_str]);
ft_plot_ortho(mritem.anatomy,'location',location_mri(1:3),'transform',inv(mritem.transform));
