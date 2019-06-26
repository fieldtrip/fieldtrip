function dartel = tbx_cfg_dartel
% Configuration file for toolbox 'Dartel Tools'
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: tbx_cfg_dartel.m 6522 2015-08-14 18:55:42Z john $

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','DARTEL')); end

% ---------------------------------------------------------------------
% matnames Parameter Files
% ---------------------------------------------------------------------
matnames         = cfg_files;
matnames.tag     = 'matnames';
matnames.name    = 'Parameter Files';
matnames.help    = {'Select ''_sn.mat'' files containing the spatial transformation and segmentation parameters. Rigidly aligned versions of the image that was segmented will be generated. The image files used by the segmentation may have moved. If they have, then (so the import can find them) ensure that they are either in the output directory, or the current working directory.'};
matnames.filter = 'mat';
matnames.ufilter = '.*seg_sn\.mat$';
matnames.num     = [1 Inf];
% ---------------------------------------------------------------------
% odir Output Directory
% ---------------------------------------------------------------------
odir         = cfg_files;
odir.tag     = 'odir';
odir.name    = 'Output Directory';
odir.help    = {'Select the directory where the resliced files should be written.'};
odir.filter = 'dir';
odir.ufilter = '.*';
odir.num     = [1 1];
% ---------------------------------------------------------------------
% bb Bounding box
% ---------------------------------------------------------------------
bb         = cfg_entry;
bb.tag     = 'bb';
bb.name    = 'Bounding box';
bb.val     = {[NaN NaN NaN
                    NaN NaN NaN]};
bb.help    = {'The bounding box (in mm) of the volume that is to be written (relative to the anterior commissure). Non-finite values will be replaced by the bounding box of the tissue probability maps used in the segmentation.'};
bb.strtype = 'e';
bb.num     = [2 3];
% ---------------------------------------------------------------------
% vox Voxel size
% ---------------------------------------------------------------------
vox         = cfg_entry;
vox.tag     = 'vox';
vox.name    = 'Voxel size';
vox.val     = {1.5};
vox.help    = {'The (isotropic) voxel sizes of the written images. A non-finite value will be replaced by the average voxel size of the tissue probability maps used by the segmentation.'};
vox.strtype = 'e';
vox.num     = [1 1];
% ---------------------------------------------------------------------
% image Image option
% ---------------------------------------------------------------------
image         = cfg_menu;
image.tag     = 'image';
image.name    = 'Image option';
image.val     = {0};
image.help    = {'A resliced version of the original image can be produced, which may have various procedures applied to it.  All options will rescale the images so that the mean of the white matter intensity is set to one. The ``skull stripped'''' versions are the images simply scaled by the sum of the grey and white matter probabilities.'};
image.labels  = {
                'Original'
                'Bias Corrected'
                'Skull-Stripped'
                'Bias Corrected and Skull-stripped'
                'None'
}';
image.values = {1 3 5 7 0};
% ---------------------------------------------------------------------
% GM Grey Matter
% ---------------------------------------------------------------------
GM         = cfg_menu;
GM.tag     = 'GM';
GM.name    = 'Grey Matter';
GM.val     = {1};
GM.help    = {'Produce a resliced version of this tissue class?'};
GM.labels  = {
             'Yes'
             'No'
}';
GM.values  = {1 0};
% ---------------------------------------------------------------------
% WM White Matter
% ---------------------------------------------------------------------
WM         = cfg_menu;
WM.tag     = 'WM';
WM.name    = 'White Matter';
WM.val     = {1};
WM.help    = {'Produce a resliced version of this tissue class?'};
WM.labels  = {
             'Yes'
             'No'
}';
WM.values  = {1 0};
% ---------------------------------------------------------------------
% CSF CSF
% ---------------------------------------------------------------------
CSF         = cfg_menu;
CSF.tag     = 'CSF';
CSF.name    = 'CSF';
CSF.val     = {0};
CSF.help    = {'Produce a resliced version of this tissue class?'};
CSF.labels  = {
              'Yes'
              'No'
}';
CSF.values  = {1 0};
% ---------------------------------------------------------------------
% initial Initial Import
% ---------------------------------------------------------------------
initial         = cfg_exbranch;
initial.tag     = 'initial';
initial.name    = 'Initial Import';
initial.val     = {matnames odir bb vox image GM WM CSF };
initial.help    = {'Images first need to be imported into a form that Dartel can work with. If the default segmentation is used (ie the Segment button), then this involves taking the results of the segmentation (*_seg_sn.mat)/* \cite{ashburner05} */, in order to have rigidly aligned tissue class images. Typically, there would be imported grey matter and white matter images, but CSF images can also be included. The subsequent Dartel alignment will then attempt to nonlinearly register these tissue class images together. If the new segmentation routine is used (from the toolbox), then this includes the option to generate ``imported'''' tissue class images. This means that a separate importing step is not needed for it.'};
initial.prog = @spm_dartel_import;
initial.vout = @vout_initial_import;
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images1         = cfg_files;
images1.tag     = 'images';
images1.name    = 'Images';
images1.help    = {'Select a set of imported images of the same type to be registered by minimising a measure of difference from the template.'};
images1.filter = 'image';
images1.ufilter = '^r.*';
images1.num     = [1 Inf];
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images         = cfg_repeat;
images.tag     = 'images';
images.name    = 'Images';
images.help    = {'Select the images to be warped together. Multiple sets of images can be simultaneously registered. For example, the first set may be a bunch of grey matter images, and the second set may be the white matter images of the same subjects.'};
images.values  = {images1 };
images.num     = [1 Inf];
% ---------------------------------------------------------------------
% template Template basename
% ---------------------------------------------------------------------
template         = cfg_entry;
template.tag     = 'template';
template.name    = 'Template basename';
template.val     = {'Template'};
template.help    = {'Enter the base for the template name.  Templates generated at each outer iteration of the procedure will be basename_1.nii, basename_2.nii etc.  If empty, then no template will be saved. Similarly, the estimated flow-fields will have the basename appended to them.'};
template.strtype = 's';
template.num     = [1 Inf];
% ---------------------------------------------------------------------
% rform Regularisation Form
% ---------------------------------------------------------------------
rform         = cfg_menu;
rform.tag     = 'rform';
rform.name    = 'Regularisation Form';
rform.val     = {0};
rform.help    = {'The registration is penalised by some ``energy'''' term.  Here, the form of this energy term is specified. Three different forms of regularisation can currently be used.'};
rform.labels = {
                'Linear Elastic Energy'
                'Membrane Energy'
                'Bending Energy'
}';
rform.values  = {0 1 2};
% ---------------------------------------------------------------------
% its Inner Iterations
% ---------------------------------------------------------------------
its         = cfg_menu;
its.tag     = 'its';
its.name    = 'Inner Iterations';
its.val     = {3};
its.help    = {'The number of Gauss-Newton iterations to be done within this outer iteration. After this, new average(s) are created, which the individual images are warped to match.'};
its.labels  = {
              '1'
              '2'
              '3'
              '4'
              '5'
              '6'
              '7'
              '8'
              '9'
              '10'
}';
its.values  = {1 2 3 4 5 6 7 8 9 10};
% ---------------------------------------------------------------------
% rparam Reg params
% ---------------------------------------------------------------------
rparam         = cfg_entry;
rparam.tag     = 'rparam';
rparam.name    = 'Reg params';
rparam.val     = {[0.1 0.01 0.001]};
rparam.help    = {
                  'For linear elasticity, the parameters are mu, lambda and id. For membrane energy, the parameters are lambda, unused and id.id is a term for penalising absolute displacements, and should therefore be small.  For bending energy, the parameters are lambda, id1 and id2, and the regularisation is by (-lambda*Laplacian + id1)^2 + id2.'
                  'Use more regularisation for the early iterations so that the deformations are smooth, and then use less for the later ones so that the details can be better matched.'
}';
rparam.strtype = 'e';
rparam.num     = [1 3];
% ---------------------------------------------------------------------
% K Time Steps
% ---------------------------------------------------------------------
K         = cfg_menu;
K.tag     = 'K';
K.name    = 'Time Steps';
K.val     = {6};
K.help    = {'The number of time points used for solving the partial differential equations.  A single time point would be equivalent to a small deformation model. Smaller values allow faster computations, but are less accurate in terms of inverse consistency and may result in the one-to-one mapping breaking down.  Earlier iteration could use fewer time points, but later ones should use about 64 (or fewer if the deformations are very smooth).'};
K.labels = {
            '1'
            '2'
            '4'
            '8'
            '16'
            '32'
            '64'
            '128'
            '256'
            '512'
}';
K.values = {0 1 2 3 4 5 6 7 8 9};
% ---------------------------------------------------------------------
% slam Smoothing Parameter
% ---------------------------------------------------------------------
slam         = cfg_menu;
slam.tag     = 'slam';
slam.name    = 'Smoothing Parameter';
slam.val     = {1};
slam.help    = {'A LogOdds parameterisation of the template is smoothed using a multi-grid scheme.  The amount of smoothing is determined by this parameter.'};
slam.labels  = {
               'None'
               '0.5'
               '1'
               '2'
               '4'
               '8'
               '16'
               '32'
}';
slam.values  = {0 0.5 1 2 4 8 16 32};
% ---------------------------------------------------------------------
% param Outer Iteration
% ---------------------------------------------------------------------
param1         = cfg_branch;
param1.tag     = 'param';
param1.name    = 'Outer Iteration';
param1.val     = {its rparam K slam };
param1.help    = {'Different parameters can be specified for each outer iteration. Each of them warps the images to the template, and then regenerates the template from the average of the warped images. Multiple outer iterations should be used for more accurate results, beginning with a more coarse registration (more regularisation) then ending with the more detailed registration (less regularisation).'};

val = cell(1,6);
param = param1;
param.val{1}.val{1} = 3; % iits
param.val{2}.val{1} = [4 2 1e-6]; % rparam
param.val{3}.val{1} = 0; % K
param.val{4}.val{1} = 16;
val{1} = param;
param.val{1}.val{1} = 3; % iits
param.val{2}.val{1} = [2 1 1e-6]; % rparam
param.val{3}.val{1} = 0; % K
param.val{4}.val{1} = 8;
val{2} = param;
param.val{1}.val{1} = 3; % iits
param.val{2}.val{1} = [1 0.5 1e-6]; % rparam
param.val{3}.val{1} = 1; % K
param.val{4}.val{1} = 4;
val{3} = param;
param.val{1}.val{1} = 3; % iits
param.val{2}.val{1} = [0.5 0.25 1e-6]; % rparam
param.val{3}.val{1} = 2; % K
param.val{4}.val{1} = 2;
val{4} = param;
param.val{1}.val{1} = 3; % iits
param.val{2}.val{1} = [0.25 0.125 1e-6]; % rparam
param.val{3}.val{1} = 4; % K
param.val{4}.val{1} = 1;
val{5} = param;
param.val{1}.val{1} = 3; % iits
param.val{2}.val{1} = [0.25 0.125 1e-6]; % rparam
param.val{3}.val{1} = 6; % K
param.val{4}.val{1} = 0.5;
val{6} = param;

% ---------------------------------------------------------------------
% param Outer Iterations
% ---------------------------------------------------------------------
param         = cfg_repeat;
param.tag     = 'param';
param.name    = 'Outer Iterations';
param.val     = val;
param.help    = {'The images are averaged, and each individual image is warped to match this average.  This is repeated a number of times.'};
param.values  = {param1};
param.num     = [1 Inf];
% ---------------------------------------------------------------------
% lmreg LM Regularisation
% ---------------------------------------------------------------------
lmreg         = cfg_entry;
lmreg.tag     = 'lmreg';
lmreg.name    = 'LM Regularisation';
lmreg.val     = {0.01};
lmreg.help    = {'Levenberg-Marquardt regularisation.  Larger values increase the the stability of the optimisation, but slow it down.  A value of zero results in a Gauss-Newton strategy, but this is not recommended as it may result in instabilities in the FMG.'};
lmreg.strtype = 'e';
lmreg.num     = [1 1];
% ---------------------------------------------------------------------
% cyc Cycles
% ---------------------------------------------------------------------
cyc         = cfg_menu;
cyc.tag     = 'cyc';
cyc.name    = 'Cycles';
cyc.val     = {3};
cyc.help    = {'Number of cycles used by the full multi-grid matrix solver. More cycles result in higher accuracy, but slow down the algorithm. See Numerical Recipes for more information on multi-grid methods.'};
cyc.labels  = {
              '1'
              '2'
              '3'
              '4'
              '5'
              '6'
              '7'
              '8'
}';
cyc.values  = {1 2 3 4 5 6 7 8};
% ---------------------------------------------------------------------
% its Iterations
% ---------------------------------------------------------------------
its         = cfg_menu;
its.tag     = 'its';
its.name    = 'Iterations';
its.val     = {3};
its.help    = {'Number of relaxation iterations performed in each multi-grid cycle. More iterations are needed if using ``bending energy'''' regularisation, because the relaxation scheme only runs very slowly. See the chapter on solving partial differential equations in Numerical Recipes for more information about relaxation methods.'};
its.labels  = {
              '1'
              '2'
              '3'
              '4'
              '5'
              '6'
              '7'
              '8'
}';
its.values   = {1 2 3 4 5 6 7 8};
% ---------------------------------------------------------------------
% optim Optimisation Settings
% ---------------------------------------------------------------------
optim         = cfg_branch;
optim.tag     = 'optim';
optim.name    = 'Optimisation Settings';
optim.val     = {lmreg cyc its };
optim.help    = {'Settings for the optimisation.  If you are unsure about them, then leave them at the default values.  Optimisation is by repeating a number of Levenberg-Marquardt iterations, in which the equations are solved using a full multi-grid (FMG) scheme. FMG and Levenberg-Marquardt are both described in Numerical Recipes (2nd edition).'};
% ---------------------------------------------------------------------
% settings Settings
% ---------------------------------------------------------------------
settings         = cfg_branch;
settings.tag     = 'settings';
settings.name    = 'Settings';
settings.val     = {template rform param optim };
settings.help    = {'Various settings for the optimisation. The default values should work reasonably well for aligning tissue class images together.'};
% ---------------------------------------------------------------------
% warp Run Dartel (create Templates)
% ---------------------------------------------------------------------
warp         = cfg_exbranch;
warp.tag     = 'warp';
warp.name    = 'Run Dartel (create Templates)';
warp.val     = {images settings };
warp.check   = @check_dartel_template;
warp.help    = {'Run the Dartel nonlinear image registration procedure. This involves iteratively matching all the selected images to a template generated from their own mean. A series of Template*.nii files are generated, which become increasingly crisp as the registration proceeds.'};
warp.prog = @spm_dartel_template;
warp.vout = @vout_dartel_template;
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images1         = cfg_files;
images1.tag     = 'images';
images1.name    = 'Images';
images1.help    = {'Select a set of imported images of the same type to be registered by minimising a measure of difference from the template.'};
images1.filter = 'image';
images1.ufilter = '^r.*';
images1.num     = [1 Inf];
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images         = cfg_repeat;
images.tag     = 'images';
images.name    = 'Images';
images.help    = {'Select the images to be warped together. Multiple sets of images can be simultaneously registered. For example, the first set may be a bunch of grey matter images, and the second set may be the white matter images of the same subjects.'};
images.values  = {images1 };
images.num     = [1 Inf];
% ---------------------------------------------------------------------
% rform Regularisation Form
% ---------------------------------------------------------------------
rform         = cfg_menu;
rform.tag     = 'rform';
rform.name    = 'Regularisation Form';
rform.val     = {0};
rform.help    = {'The registration is penalised by some ``energy'''' term.  Here, the form of this energy term is specified. Three different forms of regularisation can currently be used.'};
rform.labels  = {
                'Linear Elastic Energy'
                'Membrane Energy'
                'Bending Energy'
}';
rform.values  = {0 1 2};
% ---------------------------------------------------------------------
% its Inner Iterations
% ---------------------------------------------------------------------
its         = cfg_menu;
its.tag     = 'its';
its.name    = 'Inner Iterations';
its.val     = {3};
its.help    = {'The number of Gauss-Newton iterations to be done within this outer iteration.'};
its.labels  = {
              '1'
              '2'
              '3'
              '4'
              '5'
              '6'
              '7'
              '8'
              '9'
              '10'
}';
its.values  = {1 2 3 4 5 6 7 8 9 10};
% ---------------------------------------------------------------------
% rparam Reg params
% ---------------------------------------------------------------------
rparam         = cfg_entry;
rparam.tag     = 'rparam';
rparam.name    = 'Reg params';
rparam.val     = {[0.1 0.01 0.001]};
rparam.help    = {
                  'For linear elasticity, the parameters are mu, lambda and id. For membrane energy, the parameters are lambda, unused and id.id is a term for penalising absolute displacements, and should therefore be small.  For bending energy, the parameters are lambda, id1 and id2, and the regularisation is by (-lambda*Laplacian + id1)^2 + id2.'
                  'Use more regularisation for the early iterations so that the deformations are smooth, and then use less for the later ones so that the details can be better matched.'
}';
rparam.strtype = 'e';
rparam.num     = [1 3];
% ---------------------------------------------------------------------
% K Time Steps
% ---------------------------------------------------------------------
K         = cfg_menu;
K.tag     = 'K';
K.name    = 'Time Steps';
K.val     = {6};
K.help    = {'The number of time points used for solving the partial differential equations.  A single time point would be equivalent to a small deformation model. Smaller values allow faster computations, but are less accurate in terms of inverse consistency and may result in the one-to-one mapping breaking down.  Earlier iteration could use fewer time points, but later ones should use about 64 (or fewer if the deformations are very smooth).'};
K.labels  = {
            '1'
            '2'
            '4'
            '8'
            '16'
            '32'
            '64'
            '128'
            '256'
            '512'
}';
K.values  = {0 1 2 3 4 5 6 7 8 9};
% ---------------------------------------------------------------------
% template Template
% ---------------------------------------------------------------------
template         = cfg_files;
template.tag     = 'template';
template.name    = 'Template';
template.help    = {'Select template. Smoother templates should be used for the early iterations. Note that the template should be a 4D file, with the 4th dimension equal to the number of sets of images.'};
template.filter = 'nifti';
template.ufilter = '.*';
template.num     = [1 1];
% ---------------------------------------------------------------------
% param Outer Iteration
% ---------------------------------------------------------------------
param1         = cfg_branch;
param1.tag     = 'param';
param1.name    = 'Outer Iteration';
param1.val     = {its rparam K template };
param1.help    = {'Different parameters and templates can be specified for each outer iteration.'};

val = cell(1,6);
param = param1;
param.val{1}.val{1} = 3; % iits
param.val{2}.val{1} = [4 2 1e-6]; % rparam
param.val{3}.val{1} = 0; % K
val{1} = param;
param.val{1}.val{1} = 3; % iits
param.val{2}.val{1} = [2 1 1e-6]; % rparam
param.val{3}.val{1} = 0; % K
val{2} = param;
param.val{1}.val{1} = 3; % iits
param.val{2}.val{1} = [1 0.5 1e-6]; % rparam
param.val{3}.val{1} = 1; % K
val{3} = param;
param.val{1}.val{1} = 3; % iits
param.val{2}.val{1} = [0.5 0.25 1e-6]; % rparam
param.val{3}.val{1} = 2; % K
val{4} = param;
param.val{1}.val{1} = 3; % iits
param.val{2}.val{1} = [0.25 0.125 1e-6]; % rparam
param.val{3}.val{1} = 4; % K
val{5} = param;
param.val{1}.val{1} = 3; % iits
param.val{2}.val{1} = [0.25 0.125 1e-6]; % rparam
param.val{3}.val{1} = 6; % K
val{6} = param;

% ---------------------------------------------------------------------
% param Outer Iterations
% ---------------------------------------------------------------------
param         = cfg_repeat;
param.tag     = 'param';
param.name    = 'Outer Iterations';
param.val     = val;
param.help    = {'The images are warped to match a sequence of templates. Early iterations should ideally use smoother templates and more regularisation than later iterations.'};
param.values  = {param1};
param.num     = [1 Inf];
% ---------------------------------------------------------------------
% lmreg LM Regularisation
% ---------------------------------------------------------------------
lmreg         = cfg_entry;
lmreg.tag     = 'lmreg';
lmreg.name    = 'LM Regularisation';
lmreg.val     = {0.01};
lmreg.help    = {'Levenberg-Marquardt regularisation.  Larger values increase the the stability of the optimisation, but slow it down.  A value of zero results in a Gauss-Newton strategy, but this is not recommended as it may result in instabilities in the FMG.'};
lmreg.strtype = 'e';
lmreg.num     = [1 1];
% ---------------------------------------------------------------------
% cyc Cycles
% ---------------------------------------------------------------------
cyc         = cfg_menu;
cyc.tag     = 'cyc';
cyc.name    = 'Cycles';
cyc.val     = {3};
cyc.help    = {'Number of cycles used by the full multi-grid matrix solver. More cycles result in higher accuracy, but slow down the algorithm. See Numerical Recipes for more information on multi-grid methods.'};
cyc.labels = {
              '1'
              '2'
              '3'
              '4'
              '5'
              '6'
              '7'
              '8'
}';
cyc.values  = {1 2 3 4 5 6 7 8};
% ---------------------------------------------------------------------
% its Iterations
% ---------------------------------------------------------------------
its         = cfg_menu;
its.tag     = 'its';
its.name    = 'Iterations';
its.val     = {3};
its.help    = {'Number of relaxation iterations performed in each multi-grid cycle. More iterations are needed if using ``bending energy'''' regularisation, because the relaxation scheme only runs very slowly. See the chapter on solving partial differential equations in Numerical Recipes for more information about relaxation methods.'};
its.labels  = {
              '1'
              '2'
              '3'
              '4'
              '5'
              '6'
              '7'
              '8'
}';
its.values  = {1 2 3 4 5 6 7 8};
% ---------------------------------------------------------------------
% optim Optimisation Settings
% ---------------------------------------------------------------------
optim         = cfg_branch;
optim.tag     = 'optim';
optim.name    = 'Optimisation Settings';
optim.val     = {lmreg cyc its };
optim.help    = {'Settings for the optimisation.  If you are unsure about them, then leave them at the default values.  Optimisation is by repeating a number of Levenberg-Marquardt iterations, in which the equations are solved using a full multi-grid (FMG) scheme. FMG and Levenberg-Marquardt are both described in Numerical Recipes (2nd edition).'};
% ---------------------------------------------------------------------
% settings Settings
% ---------------------------------------------------------------------
settings         = cfg_branch;
settings.tag     = 'settings';
settings.name    = 'Settings';
settings.val     = {rform param optim };
settings.help    = {'Various settings for the optimisation. The default values should work reasonably well for aligning tissue class images together.'};
% ---------------------------------------------------------------------
% warp1 Run Dartel (existing Templates)
% ---------------------------------------------------------------------
warp1         = cfg_exbranch;
warp1.tag     = 'warp1';
warp1.name    = 'Run Dartel (existing Templates)';
warp1.val     = {images settings };
warp1.check   = @check_dartel_template;
warp1.help    = {'Run the Dartel nonlinear image registration procedure to match individual images to pre-existing template data. Start out with smooth templates, and select crisp templates for the later iterations.'};
warp1.prog = @spm_dartel_warp;
warp1.vout = @vout_dartel_warp;
% ---------------------------------------------------------------------
% flowfields Flow fields
% ---------------------------------------------------------------------
flowfields         = cfg_files;
flowfields.tag     = 'flowfields';
flowfields.name    = 'Flow fields';
flowfields.help    = {'The flow fields store the deformation information. The same fields can be used for both forward or backward deformations (or even, in principle, half way or exaggerated deformations).'};
flowfields.filter = 'nifti';
flowfields.ufilter = '^u_.*';
flowfields.num     = [1 Inf];
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images1         = cfg_files;
images1.tag     = 'images';
images1.name    = 'Images';
images1.help    = {'Select images to be warped. Note that there should be the same number of images as there are flow fields, such that each flow field warps one image.'};
images1.filter = 'nifti';
images1.ufilter = '.*';
images1.num     = [1 Inf];
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images         = cfg_repeat;
images.tag     = 'images';
images.name    = 'Images';
images.help    = {'The flow field deformations can be applied to multiple images. At this point, you are choosing how many images each flow field should be applied to.'};
images.values  = {images1 };
images.num     = [1 Inf];
% ---------------------------------------------------------------------
% See later
% ---------------------------------------------------------------------
many_subj = cfg_branch;
many_subj.tag = 'subjs';
many_subj.name = 'Many Subjects';
many_subj.val  = {flowfields,images};
many_subj.help = {[...
'Select this option if you have many subjects to spatially normalise, ',...
'but there are a small and fixed number of scans for each subject.']};
% ---------------------------------------------------------------------
% jactransf Modulation
% ---------------------------------------------------------------------
jactransf         = cfg_menu;
jactransf.tag     = 'jactransf';
jactransf.name    = 'Modulation';
jactransf.val     = {0};
jactransf.help    = {'This allows the spatially normalised images to be rescaled by the Jacobian determinants of the deformations. Note that the rescaling is only approximate for deformations generated using smaller numbers of time steps.'};
jactransf.labels  = {
                    'Pres. Concentration (No "modulation")'
                    'Pres. Amount ("Modulation")'
}';
jactransf.values  = {0 1};
% ---------------------------------------------------------------------
% K Time Steps
% ---------------------------------------------------------------------
K         = cfg_menu;
K.tag     = 'K';
K.name    = 'Time Steps';
K.val     = {6};
K.help    = {'The number of time points used for solving the partial differential equations.  Note that Jacobian determinants are not very accurate for very small numbers of time steps (less than about 16).'};
K.labels  = {
            '1'
            '2'
            '4'
            '8'
            '16'
            '32'
            '64'
            '128'
            '256'
            '512'
}';
K.values  = {0 1 2 3 4 5 6 7 8 9};
% ---------------------------------------------------------------------
% interp Interpolation
% ---------------------------------------------------------------------
interp         = cfg_menu;
interp.tag     = 'interp';
interp.name    = 'Interpolation';
interp.val     = {1};
interp.help    = {
                  ['The method by which the images are sampled when ',...
                  'being written in a different space. ',...
                  '(Note that Inf or NaN values are treated as zero, ',...
                  'rather than as missing data)']
                  '    Nearest Neighbour:'
                  '      - Fastest, but not normally recommended.'
                  '    Trilinear Interpolation:'
                  '      - OK for PET, realigned fMRI, or segmentations'
                  '    B-spline Interpolation:'
                  ['      - Better quality (but slower) interpolation',...
                  '/* \cite{thevenaz00a}*/, especially with higher ',...
                  'degree splines. Can produce values outside the ',...
                  'original range (e.g. small negative values from an ',...
                  'originally all positive image).']
}';
interp.labels  = {
                 'Nearest neighbour'
                 'Trilinear'
                 '2nd Degree B-spline'
                 '3rd Degree B-Spline '
                 '4th Degree B-Spline '
                 '5th Degree B-Spline'
                 '6th Degree B-Spline'
                 '7th Degree B-Spline'
}';
interp.values  = {0 1 2 3 4 5 6 7};
% ---------------------------------------------------------------------
% crt_warped Create Warped
% ---------------------------------------------------------------------
crt_warped         = cfg_exbranch;
crt_warped.tag     = 'crt_warped';
crt_warped.name    = 'Create Warped';
crt_warped.val     = {flowfields images jactransf K interp };
crt_warped.check   = @check_norm;
crt_warped.help    = {'This allows spatially normalised images to be generated. Note that voxel sizes and bounding boxes can not be adjusted, and that there may be strange effects due to the boundary conditions used by the warping. Also note that the warped images are not in Talairach or MNI space. The coordinate system is that of the average shape and size of the subjects to which Dartel was applied. In order to have MNI-space normalised images, then the Deformations Utility can be used to compose the individual Dartel warps, with a deformation field that matches (e.g.) the Template grey matter generated by Dartel, with one of the grey matter volumes released with SPM.'};
crt_warped.prog = @spm_dartel_norm;
crt_warped.vout = @vout_norm;
% ---------------------------------------------------------------------
% flowfields Flow fields
% ---------------------------------------------------------------------
flowfields         = cfg_files;
flowfields.tag     = 'flowfields';
flowfields.name    = 'Flow fields';
flowfields.help    = {'The flow fields store the deformation information. The same fields can be used for both forward or backward deformations (or even, in principle, half way or exaggerated deformations).'};
flowfields.filter = 'nifti';
flowfields.ufilter = '^u_.*';
flowfields.num     = [1 Inf];
% ---------------------------------------------------------------------
% K Time Steps
% ---------------------------------------------------------------------
K         = cfg_menu;
K.tag     = 'K';
K.name    = 'Time Steps';
K.val     = {6};
K.help    = {'The number of time points used for solving the partial differential equations.  Note that Jacobian determinants are not very accurate for very small numbers of time steps (less than about 16).'};
K.labels  = {
            '1'
            '2'
            '4'
            '8'
            '16'
            '32'
            '64'
            '128'
            '256'
            '512'
}';
K.values = {0 1 2 3 4 5 6 7 8 9};
% ---------------------------------------------------------------------
% jacdet Jacobian determinants
% ---------------------------------------------------------------------
jacdet         = cfg_exbranch;
jacdet.tag     = 'jacdet';
jacdet.name    = 'Jacobian determinants';
jacdet.val     = {flowfields K };
jacdet.help    = {'Create Jacobian determinant fields from flowfields.'};
jacdet.prog = @spm_dartel_jacobian;
jacdet.vout = @vout_jacdet;
% ---------------------------------------------------------------------
% flowfields Flow fields
% ---------------------------------------------------------------------
flowfields         = cfg_files;
flowfields.tag     = 'flowfields';
flowfields.name    = 'Flow fields';
flowfields.help    = {'The flow fields store the deformation information. The same fields can be used for both forward or backward deformations (or even, in principle, half way or exaggerated deformations).'};
flowfields.filter = 'nifti';
flowfields.ufilter = '^u_.*';
flowfields.num     = [1 Inf];
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images         = cfg_files;
images.tag     = 'images';
images.name    = 'Images';
images.help    = {'Select the image(s) to be inverse normalised.  These should be in alignment with the template image of the warping procedure (Run Dartel).'};
images.filter = 'nifti';
images.ufilter = '.*';
images.num     = [1 Inf];
% ---------------------------------------------------------------------
% K Time Steps
% ---------------------------------------------------------------------
K         = cfg_menu;
K.tag     = 'K';
K.name    = 'Time Steps';
K.val     = {6};
K.help    = {'The number of time points used for solving the partial differential equations.  Note that Jacobian determinants are not very accurate for very small numbers of time steps (less than about 16).'};
K.labels  = {
            '1'
            '2'
            '4'
            '8'
            '16'
            '32'
            '64'
            '128'
            '256'
            '512'
}';
K.values = {0 1 2 3 4 5 6 7 8 9};
% ---------------------------------------------------------------------
% crt_iwarped Create Inverse Warped
% ---------------------------------------------------------------------
crt_iwarped         = cfg_exbranch;
crt_iwarped.tag     = 'crt_iwarped';
crt_iwarped.name    = 'Create Inverse Warped';
crt_iwarped.val     = {flowfields images K interp };
crt_iwarped.help    = {'Create inverse normalised versions of some image(s). The image that is inverse-normalised should be in alignment with the template (generated during the warping procedure). Note that the results have the same dimensions as the ``flow fields'''', but are mapped to the original images via the affine transformations in their headers.'};
crt_iwarped.prog = @spm_dartel_invnorm;
crt_iwarped.vout = @vout_invnorm;
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
flowfield         = cfg_files;
flowfield.tag     =  'flowfield';
flowfield.name    = 'Flow Field';
flowfield.filter  = 'nifti';
flowfield.ufilter = '^u_.*\.nii$';
flowfield.num     = [1 1];
flowfield.help    = {'Dartel flow field for this subject.'};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
images        = cfg_files;
images.tag    = 'images';
images.name   = 'Images';
images.filter = 'nifti';
images.num    = [1 Inf];
images.help   = {'Images for this subject to spatially normalise.'};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
subj      = cfg_branch;
subj.tag  = 'subj';
subj.name = 'Subject';
subj.val  = {flowfield,images};
subj.help = {'Subject to be spatially normalized.'};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
few_subj        = cfg_repeat;
few_subj.tag    = 'few_subj';
few_subj.name   = 'Few Subjects';
few_subj.values = {subj};
few_subj.help   = {[...
'Select this option if there are only a few subjects, each with many or ' ...
'a variable number of scans each. You will then need to specify a series of subjects, and the flow field and images of each of them.']};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
way        = cfg_choice;
way.tag    = 'data';
way.name   = 'Select according to';
way.values = {few_subj,many_subj};
way.help   = {...
['You may wish to spatially normalise only a few subjects, '...
 'but have many scans per subject (eg for fMRI), '...
 'or you may have lots of subjects, but with a small and fixed number '...
 'of scans for each of them (eg for VBM).  The idea is to chose the way of '...
 'selecting files that is easier.']};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
template        = cfg_files;
template.tag    = 'template';
template.name   = 'Dartel Template';
template.filter = 'nifti';
template.num    = [0 1];
template.help   = {...
['Select the final Template file generated by Dartel. This will be affine '...
 'registered with a TPM file, such that the resulting spatially normalised '...
 'images are closer aligned to MNI space. Leave empty if you do not wish to '...
 'incorporate a transform to MNI space '...
 '(ie just click ``done'' on the file selector, without selecting any images).']};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'Gaussian FWHM';
fwhm.val     = {[8 8 8]};
fwhm.strtype = 'e';
fwhm.num     = [1 3];
fwhm.help    = {'Specify the full-width at half maximum (FWHM) of the Gaussian blurring kernel in mm. Three values should be entered, denoting the FWHM in the x, y and z directions. Note that you can also specify [0 0 0], but any ``modulated'' data will show aliasing (see eg Wikipedia), which occurs because of the way the warped images are generated.'};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
preserve         = cfg_menu;
preserve.tag     = 'preserve';
preserve.name    = 'Preserve';
preserve.help    = {
'Preserve Concentrations (no "modulation"): Smoothed spatially normalised images (sw*) represent weighted averages of the signal under the smoothing kernel, approximately preserving the intensities of the original images. This option is currently suggested for eg fMRI.'
''
'Preserve Amount ("modulation"): Smoothed and spatially normalised images preserve the total amount of signal from each region in the images (smw*). Areas that are expanded during warping are correspondingly reduced in intensity. This option is suggested for VBM.'
}';
preserve.labels = {
                   'Preserve Concentrations'
                   'Preserve Amount'
}';
preserve.values = {0 1};
preserve.val    = {0};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
vox          = cfg_entry;
vox.tag      = 'vox';
vox.name     = 'Voxel sizes';
vox.num      = [1 3];
vox.strtype  = 'e';
vox.val      = {[NaN NaN NaN]};
vox.help     = {[...
'Specify the voxel sizes of the deformation field to be produced. ',...
'Non-finite values will default to the voxel sizes of the template image',...
'that was originally used to estimate the deformation.']};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
bb           = cfg_entry;
bb.tag       = 'bb';
bb.name      = 'Bounding box';
bb.strtype   = 'e';
bb.num       = [2 3];
bb.val       = {[NaN NaN NaN; NaN NaN NaN]};
bb.help      = {[...
'Specify the bounding box of the deformation field to be produced. ',...
'Non-finite values will default to the bounding box of the template image',...
'that was originally used to estimate the deformation.']};
% ---------------------------------------------------------------------
%
% ---------------------------------------------------------------------
nrm       = cfg_exbranch;
nrm.tag   = 'mni_norm';
nrm.name  = 'Normalise to MNI Space';
nrm.val   = {template,way,vox,bb,preserve,fwhm};
nrm.prog  = @spm_dartel_norm_fun;
nrm.vout  = @vout_norm_fun;
nrm.check = @check_norm_fun;
nrm.help  = {[...
'Normally, Dartel generates warped images that align with the average-shaped template. ',...
'This routine includes an initial affine regisration of the template (the final one ',...
'generated by Dartel), with the TPM data released with SPM.'],[...
'``Smoothed'''' (blurred) spatially normalised images are generated in such a ',...
'way that the original signal is preserved. Normalised images are ',...
'generated by a ``pushing'''' rather than a ``pulling'''' (the usual) procedure. ',...
'Note that a procedure related to trilinear interpolation is used, and no masking is done.  It ',...
'is therefore recommended that the images are realigned and resliced ',...
'before they are spatially normalised, in order to benefit from motion correction using higher order interpolation.  Alternatively, contrast images ',...
'generated from unsmoothed native-space fMRI/PET data can be spatially ',...
'normalised for a 2nd level analysis.'],[...
'Two ``preserve'''' options are provided.  One of them should do the ',...
'equavalent of generating smoothed ``modulated'''' spatially normalised ',...
'images.  The other does the equivalent of smoothing the modulated ',...
'normalised fMRI/PET, and dividing by the smoothed Jacobian determinants.']};

% ---------------------------------------------------------------------
% template Dartel Template
% ---------------------------------------------------------------------
template        = cfg_files;
template.tag    = 'template';
template.name   = 'Dartel Template';
template.filter = 'nifti';
template.num    = [0 1];
template.help   = {'Select the final Template file generated by Dartel.'};
% ---------------------------------------------------------------------
% popnorm Population to ICBM Registration
% ---------------------------------------------------------------------
popnorm         = cfg_exbranch;
popnorm.tag     = 'popnorm';
popnorm.name    = 'Population to ICBM Registration';
popnorm.val     = {template};
popnorm.help    = {'Estimate the mapping from the population average to ICBM space.  This is achieved by aligning the population average with grey and white matter tissue probability maps generated by Vladimir Fonov of the McConnell Brain Imaging Centre, Montreal Neurological Institute, McGill University /* \cite{fonov2011unbiased} */.'};
popnorm.prog    = @spm_norm_population;
popnorm.vout    = @vout_popnorm;

% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images1         = cfg_files;
images1.tag     = 'images';
images1.name    = 'Images';
images1.help    = {'Select tissue class images (one per subject).'};
images1.filter = 'nifti';
images1.ufilter = '^r.*c';
images1.num     = [1 Inf];
% ---------------------------------------------------------------------
% images Images
% ---------------------------------------------------------------------
images         = cfg_repeat;
images.tag     = 'images';
images.name    = 'Images';
images.help    = {'Multiple sets of images are used here. For example, the first set may be a bunch of grey matter images, and the second set may be the white matter images of the same subjects.  The number of sets of images must be the same as was used to generate the template.'};
images.values  = {images1 };
images.num     = [1 Inf];
% ---------------------------------------------------------------------
% flowfields Flow fields
% ---------------------------------------------------------------------
flowfields         = cfg_files;
flowfields.tag     = 'flowfields';
flowfields.name    = 'Flow fields';
flowfields.help    = {'Select the flow fields for each subject.'};
flowfields.filter = 'nifti';
flowfields.ufilter = '^u_.*';
flowfields.num     = [1 Inf];
% ---------------------------------------------------------------------
% template Template
% ---------------------------------------------------------------------
template         = cfg_files;
template.tag     = 'template';
template.name    = 'Template';
template.help    = {'Residual differences are computed between the warped images and template.'};
template.filter = 'nifti';
template.ufilter = '^Template.*';
template.num     = [0 1];
% ---------------------------------------------------------------------
% K Time Steps
% ---------------------------------------------------------------------
K         = cfg_menu;
K.tag     = 'K';
K.name    = 'Time Steps';
K.val     = {6};
K.help    = {'The number of time points used for solving the partial differential equations.  Note that Jacobian determinants are not very accurate for very small numbers of time steps (less than about 16).'};
K.labels  = {
            '1'
            '2'
            '4'
            '8'
            '16'
            '32'
            '64'
            '128'
            '256'
            '512'
}';
K.values  = {0 1 2 3 4 5 6 7 8 9};
% ---------------------------------------------------------------------
% fwhm Smoothing
% ---------------------------------------------------------------------
fwhm         = cfg_menu;
fwhm.tag     = 'fwhm';
fwhm.name    = 'Smoothing';
fwhm.val     = {4};
fwhm.help    = {'The residuals can be smoothed with a Gaussian to reduce dimensionality. More smoothing is recommended if there are fewer training images.'};
fwhm.labels  = {
               'None'
               ' 2mm'
               ' 4mm'
               ' 6mm'
               ' 8mm'
               '10mm'
               '12mm'
               '14mm'
               '16mm'
}';
fwhm.values  = {0 2 4 6 8 10 12 14 16};
% ---------------------------------------------------------------------
% resids Generate Residuals
% ---------------------------------------------------------------------
resids         = cfg_exbranch;
resids.tag     = 'resids';
resids.name    = 'Generate Residuals';
resids.val     = {images flowfields template K fwhm };
resids.check   = @check_resids;
resids.help    = {'Generate residual images in a form suitable for computing a Fisher kernel. In principle, a Gaussian Process model can be used to determine the optimal (positive) linear combination of kernel matrices.  The idea would be to combine the kernel from the residuals, with a kernel derived from the flow-fields. Such a combined kernel should then encode more relevant information than the individual kernels alone.'};
resids.prog = @spm_dartel_resids;
resids.vout = @vout_resids;
% ---------------------------------------------------------------------
% images Data
% ---------------------------------------------------------------------
images         = cfg_files;
images.tag     = 'images';
images.name    = 'Data';
images.help    = {'Select images to generate dot-products from.'};
images.filter = 'nifti';
images.ufilter = '.*';
images.num     = [1 Inf];
% ---------------------------------------------------------------------
% weight Weighting image
% ---------------------------------------------------------------------
weight         = cfg_files;
weight.tag     = 'weight';
weight.name    = 'Weighting image';
weight.val = {{}};
weight.help    = {'The kernel can be generated so that some voxels contribute to the similarity measures more than others.  This is achieved by supplying a weighting image, which each of the component images are multiplied before the dot-products are computed. This image needs to have the same dimensions as the component images, but orientation information (encoded by matrices in the headers) is ignored. If left empty, then all voxels are weighted equally.'};
weight.filter = 'image';
weight.ufilter = '.*';
weight.num     = [0 1];
% ---------------------------------------------------------------------
% dotprod Dot-product Filename
% ---------------------------------------------------------------------
dotprod         = cfg_entry;
dotprod.tag     = 'dotprod';
dotprod.name    = 'Dot-product Filename';
dotprod.help    = {'Enter a filename for results (it will be prefixed by ``dp_'''' and saved in the current directory).'};
dotprod.strtype = 's';
dotprod.num     = [1 Inf];
% ---------------------------------------------------------------------
% reskern Kernel from Resids
% ---------------------------------------------------------------------
reskern         = cfg_exbranch;
reskern.tag     = 'reskern';
reskern.name    = 'Kernel from Images';
reskern.val     = {images weight dotprod };
reskern.help    = {'Generate a kernel matrix from images. In principle, this same function could be used for generating kernels from any image data (e.g. ``modulated'''' grey matter). If there is prior knowledge about some region providing more predictive information (e.g. the hippocampi for AD), then it is possible to weight the generation of the kernel accordingly. The matrix of dot-products is saved in a variable ``Phi'''', which can be loaded from the dp_*.mat file. The ``kernel trick'''' can be used to convert these dot-products into distance measures for e.g. radial basis-function approaches.'};
reskern.prog = @spm_dartel_dotprods;
% ---------------------------------------------------------------------
% flowfields Flow fields
% ---------------------------------------------------------------------
flowfields         = cfg_files;
flowfields.tag     = 'flowfields';
flowfields.name    = 'Flow fields';
flowfields.help    = {'Select the flow fields for each subject.'};
flowfields.filter = 'nifti';
flowfields.ufilter = '^u_.*';
flowfields.num     = [1 Inf];
% ---------------------------------------------------------------------
% rform Regularisation Form
% ---------------------------------------------------------------------
rform         = cfg_menu;
rform.tag     = 'rform';
rform.name    = 'Regularisation Form';
rform.val     = {0};
rform.help    = {'The registration is penalised by some ``energy'''' term.  Here, the form of this energy term is specified. Three different forms of regularisation can currently be used.'};
rform.labels  = {
                'Linear Elastic Energy'
                'Membrane Energy'
                'Bending Energy'
}';
rform.values  = {0 1 2};
% ---------------------------------------------------------------------
% rparam Reg params
% ---------------------------------------------------------------------
rparam         = cfg_entry;
rparam.tag     = 'rparam';
rparam.name    = 'Reg params';
rparam.val     = {[0.25 0.125 1e-06]};
rparam.help    = {'For linear elasticity, the parameters are `mu'', `lambda'' and `id''. For membrane and bending energy, the parameters are `lambda'', unused and `id''. The term `id'' is for penalising absolute displacements, and should therefore be small.'};
rparam.strtype = 'e';
rparam.num     = [1 3];
% ---------------------------------------------------------------------
% dotprod Dot-product Filename
% ---------------------------------------------------------------------
dotprod         = cfg_entry;
dotprod.tag     = 'dotprod';
dotprod.name    = 'Dot-product Filename';
dotprod.help    = {'Enter a filename for results (it will be prefixed by ``dp_'''' and saved in the current directory.'};
dotprod.strtype = 's';
dotprod.num     = [1 Inf];
% ---------------------------------------------------------------------
% flokern Kernel from Flows
% ---------------------------------------------------------------------
flokern         = cfg_exbranch;
flokern.tag     = 'flokern';
flokern.name    = 'Kernel from Flows';
flokern.val     = {flowfields rform rparam dotprod };
flokern.help    = {'Generate a kernel from flow fields. The dot-products are saved in a variable ``Phi'''' in the resulting dp_*.mat file.'};
flokern.prog = @spm_dartel_kernel;
% ---------------------------------------------------------------------
% kernfun Kernel Utilities
% ---------------------------------------------------------------------
kernfun         = cfg_choice;
kernfun.tag     = 'kernfun';
kernfun.name    = 'Kernel Utilities';
kernfun.help    = {
                   'Dartel can be used for generating matrices of dot-products for various kernel pattern-recognition procedures.'
                   'The idea of applying pattern-recognition procedures is to obtain a multi-variate characterisation of the anatomical differences among groups of subjects. These characterisations can then be used to separate (eg) healthy individuals from particular patient populations. There is still a great deal of methodological work to be done, so the types of kernel that can be generated here are unlikely to be the definitive ways of proceeding.  They are only just a few ideas that may be worth trying out. The idea is simply to attempt a vaguely principled way to combine generative models with discriminative models (see the ``Pattern Recognition and Machine Learning'''' book by Chris Bishop for more ideas). Better ways (higher predictive accuracy) will eventually emerge.'
                   'Various pattern recognition algorithms are available freely over the Internet. Possible approaches include Support-Vector Machines, Relevance-Vector machines and Gaussian Process Models. Gaussian Process Models probably give the most accurate probabilistic predictions, and allow kernels generated from different pieces of data to be most easily combined.'
}';
kernfun.values  = {reskern flokern};
% ---------------------------------------------------------------------
% dartel Dartel Tools
% ---------------------------------------------------------------------
dartel         = cfg_choice;
dartel.tag     = 'dartel';
dartel.name    = 'Dartel Tools';
dartel.help    = {
                  'This toolbox is based around the ``A Fast Diffeomorphic Registration Algorithm'''' paper/* \cite{ashburner07} */. The idea is to register images by computing a ``flow field'''', which can then be ``exponentiated'''' to generate both forward and backward deformations. Currently, the software only works with images that have isotropic voxels, identical dimensions and which are in approximate alignment with each other. One of the reasons for this is that the approach assumes circulant boundary conditions, which makes modelling global rotations impossible. Another reason why the images should be approximately aligned is because there are interactions among the transformations that are minimised by beginning with images that are already almost in register. This problem could be alleviated by a time varying flow field, but this is currently computationally impractical.'
                  'Because of these limitations, images should first be imported. This involves taking the ``*_seg_sn.mat'''' files produced by the segmentation code of SPM5, and writing out rigidly transformed versions of the tissue class images, such that they are in as close alignment as possible with the tissue probability maps. Rigidly transformed original images can also be generated, with the option to have skull-stripped versions.'
                  'The next step is the registration itself.  This can involve matching single images together, or it can involve the simultaneous registration of e.g. GM with GM, WM with WM and 1-(GM+WM) with 1-(GM+WM) (when needed, the 1-(GM+WM) class is generated implicitly, so there is no need to include this class yourself). This procedure begins by creating a mean of all the images, which is used as an initial template. Deformations from this template to each of the individual images are computed, and the template is then re-generated by applying the inverses of the deformations to the images and averaging. This procedure is repeated a number of times.'
                  'Finally, warped versions of the images (or other images that are in alignment with them) can be generated. '
                  ''
                  'This toolbox is not yet seamlessly integrated into the SPM package. Eventually, the plan is to use many of the ideas here as the default strategy for spatial normalisation. The toolbox may change with future updates.  There will also be a number of other (as yet unspecified) extensions, which may include a variable velocity version (related to LDDMM). Note that the Fast Diffeomorphism paper only describes a sum of squares objective function. The multinomial objective function is an extension, based on a more appropriate model for aligning binary data to a template.'
}';
dartel.values  = {initial warp warp1 nrm crt_warped jacdet crt_iwarped popnorm kernfun };
%dartel.num     = [0 Inf];

%_______________________________________________________________________
%
%_______________________________________________________________________

function dep = vout_initial_import(job)
cls = {'GM', 'WM', 'CSF'};
kk = 1;
for k=1:3,
    if isnumeric(job.(cls{k})) && job.(cls{k})
        dep(kk)            = cfg_dep;
        dep(kk).sname      = sprintf('Imported Tissue (%s)', cls{k});
        dep(kk).src_output = substruct('.','cfiles','()',{':',k});
        dep(kk).tgt_spec   = cfg_findspec({{'filter','nifti'}});
        kk = kk + 1;
    end
end
if isnumeric(job.image) && job.image
    dep(kk)            = cfg_dep;
    dep(kk).sname      = sprintf('Resliced Original Images');
    dep(kk).src_output = substruct('.','files');
    dep(kk).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
%_______________________________________________________________________

%_______________________________________________________________________
function chk = check_dartel_template(job)
n1 = numel(job.images);
n2 = numel(job.images{1});
chk = '';
for i=1:n1,
    if numel(job.images{i}) ~= n2,
        chk = 'Incompatible number of images';
        break;
    end;
end;
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_dartel_template(job)

if isa(job.settings.template,'cfg_dep') || ~ ...
        isempty(deblank(job.settings.template))
    for it=0:numel(job.settings.param),
        tdep(it+1)            = cfg_dep;
        tdep(it+1).sname      = sprintf('Template (Iteration %d)', it);
        tdep(it+1).src_output = substruct('.','template','()',{it+1});
        tdep(it+1).tgt_spec   = cfg_findspec({{'filter','nifti'}});
    end
else
    tdep = cfg_dep;
    tdep = tdep(false);
end
fdep            = cfg_dep;
fdep.sname      = 'Flow Fields';
fdep.src_output = substruct('.','files','()',{':'});
fdep.tgt_spec   = cfg_findspec({{'filter','nifti'}});
dep = [tdep fdep];
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_dartel_warp(job)
dep            = cfg_dep;
dep.sname      = 'Flow Fields';
dep.src_output = substruct('.','files','()',{':'});
dep.tgt_spec   = cfg_findspec({{'filter','nifti'}});
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_popnorm(job)
dep            = cfg_dep;
dep.sname      = 'Pop2MNI';
dep.src_output = substruct('.','files','()',{':'});
dep.tgt_spec   = cfg_findspec({{'filter','nifti'}});
%_______________________________________________________________________

%_______________________________________________________________________
function chk = check_norm(job)
chk = '';
PU = job.flowfields;
PI = job.images;
n1 = numel(PU);
for i=1:numel(PI),
    if numel(PI{i}) ~= n1,
        chk = 'Incompatible number of images';
        break;
    end
end
%_______________________________________________________________________

%_______________________________________________________________________
function chk = check_norm_fun(job)
chk = '';
if isfield(job.data,'subjs')
    PU = job.data.subjs.flowfields;
    PI = job.data.subjs.images;
    n1 = numel(PU);
    for i=1:numel(PI),
        if numel(PI{i}) ~= n1,
            chk = 'Incompatible number of images';
            break;
        end
    end
end
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_norm(job)
if job.jactransf,
    sname = 'Warped Images - Jacobian Transformed';
else
    sname = 'Warped Images';
end
PU    = job.flowfields;
PI    = job.images;
for m=1:numel(PI),
    dep(m)            = cfg_dep;
    dep(m).sname      = sprintf('%s (Image %d)',sname,m);
    dep(m).src_output = substruct('.','files','()',{':',m});
    dep(m).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_invnorm(job)
PU    = job.flowfields;
PI    = job.images;

for m=1:numel(PI),
    dep(m)            = cfg_dep;
    dep(m).sname      = sprintf('Inverse Warped Images (Image %d)',m);
    dep(m).src_output = substruct('.','files','()',{':',m});
    dep(m).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
n = numel(PI);
for m=1:numel(PU),
    dep(m+n)            = cfg_dep;
    dep(m+n).sname      = sprintf('Inverse Warped Images (Deformation %d)',m);
    dep(m+n).src_output = substruct('.','files','()',{m,':'});
    dep(m+n).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_jacdet(job)
dep            = cfg_dep;
dep.sname      = 'Jacobian Determinant Fields';
dep.src_output = substruct('.','files','()',{':'});
dep.tgt_spec   = cfg_findspec({{'filter','nifti'}});
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_norm_fun(job)
if job.preserve,
    sname = 'MNI Smo. Warped - Amount';
else
    sname = 'MNI Smo. Warped - Concentrations';
end

dep = cfg_dep;

if isfield(job.data,'subj')
    for m=1:numel(job.data.subj)
        dep(m)            = cfg_dep;
        dep(m).sname      = sprintf('%s (Deformation %d)',sname,m);
        dep(m).src_output = substruct('{}',{m},'()',{':'});
        dep(m).tgt_spec   = cfg_findspec({{'filter','nifti'}});
    end
end

if isfield(job.data,'subjs')
    for m=1:numel(job.data.subjs.images),
        dep(m)            = cfg_dep;
        dep(m).sname      = sprintf('%s (Image %d)',sname,m);
        dep(m).src_output = substruct('()',{':',m});
        dep(m).tgt_spec   = cfg_findspec({{'filter','nifti'}});
    end
end
%_______________________________________________________________________

%_______________________________________________________________________
function chk = check_resids(job)
chk = '';
PU = job.flowfields;
PI = job.images;
n1 = numel(PU);
for i=1:numel(PI),
    if numel(PI{i}) ~= n1,
        chk = 'Incompatible number of images';
        break;
    end
end
%_______________________________________________________________________

%_______________________________________________________________________
function dep = vout_resids(job)
dep = cfg_dep;
dep.sname      = 'Residual Files';
dep.src_output = substruct('.','files','()',{':'});
dep.tgt_spec   = cfg_findspec({{'filter','nifti'}});
%_______________________________________________________________________




