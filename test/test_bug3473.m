function test_bug3453
%filename = '/home/common/matlab/fieldtrip/data/ftp/tutorial/epilepsy/case1/ctf_data/case1.pos';
filename = 'D:\Development\Data\case1\ctf_data/case1.pos';
%% Individualizing a template volume conduction to on the basis of surface information
%
%% # Introduction
%
% This tutorial demonstrates how to construct an individualized template
% volume conduction model on the basis of a surface scan.
%
% This tutorial does not cover how to do the source estimation itself.
%
%% ##  Background
%
% The quality of EEG source estimates depends on accurate volume conduction
% models and sensor positions. The volume conduction model comprises a description of the geometry, of the conductivities and of a computational approach for solving Poisson's equations.
%
% The current golden standard is to measure the head geometry with an MRI
% and on the basis of that to create a volume conduction model. However,
% this data is not always available. Therefore, the usage of surface scan's
% can come in hand. In this tutorial we will describe how to use surface
% scan to create an individualized template volume conduction model.
%
%% ## Download
%
% For this the tutorial we will use
% [this](ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/epilepsy)
% dataset. More information on this dataset can be found
% [here](/tutorial/epilepsy/). 
% Also we need an external toolbox which can be downloaded [here](https://sites.google.com/site/myronenko/research/cpd)
% 
%
%
%% # Loading and coregistering data
%
% Before starting with FieldTrip, it is important that you set up your
% [MATLAB path](/faq/should_i_add_fieldtrip_with_all_subdirectories_to_my_matlab_path) properly.
%
% Then you can load the data head shape measured with the Polhemus and a
% template volume conduction model. Also we will convert the units into mm
% by now.
polhemus = ft_read_headshape(filename);
polhemus = ft_convert_units(polhemus,'mm');

template = ft_read_vol('standard_bem.mat');
template = ft_convert_units(template,'mm');
%% # Coregistration
%
% In the next step we coregister both meshes with each other.
cfg = [];
cfg.template.headshape      = polhemus;
cfg.checksize               = inf;
cfg.individual.headmodel    = template;
cfg                         = ft_interactiverealign(cfg);
template                    = ft_transform_geometry(cfg.m,template);
%% # Create surface meshes with shared features
%
% For creating the individualized mesh it is important that the head
% surface of template only contains features that are also in the head
% surface measurement of the Polhemus. Therefore, we use ft_defacemesh to remove the undesired features.

defaced_template                = template;
cfg                             = [];
defaced_template.bnd(1).unit    = 'mm';
defaced                         =  ft_defacemesh(cfg,defaced_template.bnd(1));

defaced_template.bnd(1).pos = defaced.pos;
defaced_template.bnd(1).tri = defaced.tri;
%% # Fitting template to Polhemus measurement
% 
% We will now use the surface information of the template model and the
% Polhemus measurement to create an individualised version mesh of the template mesh.
cfg             = [];
cfg.headshape   = polhemus;
cfg.template    = defaced_template.bnd(1);
cfg.method      = 'fittemplate';
fitted          = ft_prepare_mesh(cfg, template.bnd);

%% # Creating volume conduction model
%
% Finally we will create a volume conduction model.
cfg = [];
cfg.method = 'bemcp';
headmodel_bem = ft_prepare_headmodel(cfg, fitted);

%% # Summary and further reading
%
% In this tutorial we learned how to individualize a volume conduction
% model.
%
% For further reading suggest to read to about
% [BEM](/tutorial/headmodel_eeg_bem) or [FEM](/tutorial/headmodel_eeg_fem).
%
% This tutorial was last tested on 02-04-2019 by Simon Homölle on Windows 10, Matlab 2018a.

