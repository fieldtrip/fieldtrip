function mars = tbx_cfg_mars
% MATLABBATCH Configuration file for toolbox
% 'MARS (Morphologically and Anatomically accuRate Segmentation)'
%
% Adapted from tbx_cfg_preproc8 (the batch config file for toolbox
% 'New Segment') by
%
% John Ashburner
% $Id: tbx_cfg_preproc8.m 4337 2011-05-31 16:59:44Z john $
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
%
% Yu (Andy) Huang
% $Id: tbx_cfg_mars.m 2015-07-16 andy$
% Neural Engineering Lab, Dept. of Biomedical Engineering, City College of New York
% yhuang16@citymail.cuny.edu

if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','mars')); end

% ---------------------------------------------------------------------
% vols Volumes
% ---------------------------------------------------------------------
vols         = cfg_files;
vols.tag     = 'vols';
vols.name    = 'Volumes';
vols.help    = {'Select scans from this channel for processing. If multiple channels are used (eg T1 & T2), then the same order of subjects must be specified for each channel and they must be in register (same position, size, voxel dims etc..).'};
vols.filter = 'image';
vols.ufilter = '.*';
vols.num     = [1 Inf];
% ---------------------------------------------------------------------
% biasreg Bias regularisation
% ---------------------------------------------------------------------
biasreg         = cfg_menu;
biasreg.tag     = 'biasreg';
biasreg.name    = 'Bias regularisation';
biasreg.help    = {
                   'MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.'
                   ''
                   'An important issue relates to the distinction between intensity variations that arise because of bias artifact due to the physics of MR scanning, and those that arise due to different tissue properties.  The objective is to model the latter by different tissue classes, while modelling the former with a bias field. We know a priori that intensity variations due to MR physics tend to be spatially smooth, whereas those due to different tissue types tend to contain more high frequency information. A more accurate estimate of a bias field can be obtained by including prior knowledge about the distribution of the fields likely to be encountered by the correction algorithm. For example, if it is known that there is little or no intensity non-uniformity, then it would be wise to penalise large values for the intensity non-uniformity parameters. This regularisation can be placed within a Bayesian context, whereby the penalty incurred is the negative logarithm of a prior probability for any particular pattern of non-uniformity.'
                   'Knowing what works best should be a matter of empirical exploration.  For example, if your data has very little intensity non-uniformity artifact, then the bias regularisation should be increased.  This effectively tells the algorithm that there is very little bias in your data, so it does not try to model it.'
                   }';
biasreg.labels = {
                  'no regularisation (0)'
                  'extremely light regularisation (0.00001)'
                  'very light regularisation (0.0001)'
                  'light regularisation (0.001)'
                  'medium regularisation (0.01)'
                  'heavy regularisation (0.1)'
                  'very heavy regularisation (1)'
                  'extremely heavy regularisation (10)'
                  }';
biasreg.values = {
                  0
                  1e-05
                  0.0001
                  0.001
                  0.01
                  0.1
                  1
                  10
                  }';
biasreg.val    = {0.0001};
% ---------------------------------------------------------------------
% biasfwhm Bias FWHM
% ---------------------------------------------------------------------
biasfwhm         = cfg_menu;
biasfwhm.tag     = 'biasfwhm';
biasfwhm.name    = 'Bias FWHM';
biasfwhm.help    = {'FWHM of Gaussian smoothness of bias. If your intensity non-uniformity is very smooth, then choose a large FWHM. This will prevent the algorithm from trying to model out intensity variation due to different tissue types. The model for intensity non-uniformity is one of i.i.d. Gaussian noise that has been smoothed by some amount, before taking the exponential. Note also that smoother bias fields need fewer parameters to describe them. This means that the algorithm is faster for smoother intensity non-uniformities.'};
biasfwhm.labels = {
                   '30mm cutoff'
                   '40mm cutoff'
                   '50mm cutoff'
                   '60mm cutoff'
                   '70mm cutoff'
                   '80mm cutoff'
                   '90mm cutoff'
                   '100mm cutoff'
                   '110mm cutoff'
                   '120mm cutoff'
                   '130mm cutoff'
                   '140mm cutoff'
                   '150mm cutoff'
                   'No correction'
                   }';
biasfwhm.values = {
                   30
                   40
                   50
                   60
                   70
                   80
                   90
                   100
                   110
                   120
                   130
                   140
                   150
                   Inf
                   }';
biasfwhm.val    = {60};
% ---------------------------------------------------------------------
% write Save Bias Corrected
% ---------------------------------------------------------------------
write         = cfg_menu;
write.tag     = 'write';
write.name    = 'Save Bias Corrected';
write.help    = {'This is the option to save a bias corrected version of your images from this channel, or/and the estimated bias field. MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.  The bias corrected version should have more uniform intensities within the different types of tissues.'};
write.labels = {
                'Save Nothing'
                'Save Bias Corrected'
                'Save Bias Field'
                'Save Field and Corrected'
                }';
write.values = {
                [0 0]
                [0 1]
                [1 0]
                [1 1]
                }';
write.val    = {[0 0]};
% ---------------------------------------------------------------------
% channel Channel
% ---------------------------------------------------------------------
channel         = cfg_branch;
channel.tag     = 'channel';
channel.name    = 'Channel';
channel.val     = {vols biasreg biasfwhm write };
channel.help    = {'Specify a channel for processing. If multiple channels are used (eg PD & T2), then the same order of subjects must be specified for each channel and they must be in register (same position, size, voxel dims etc..). The different channels can be treated differently in terms of inhomogeneity correction etc. You may wish to correct some channels and save the corrected images, whereas you may wish not to do this for other channels.'};
% ---------------------------------------------------------------------
% data Data
% ---------------------------------------------------------------------
data         = cfg_repeat;
data.tag     = 'data';
data.name    = 'Data';
data.val     = {channel };
data.help    = {'Specify the number of different channels (for multi-spectral classification). If you have scans of different contrasts for each of the subjects, then it is possible to combine the information from them in order to improve the segmentation accuracy. Note that only the first channel of data is used for the initial affine registration with the tissue probability maps.'};
data.values  = {channel };
data.num     = [1 Inf];
% ---------------------------------------------------------------------
% tpm Tissue probability map
% ---------------------------------------------------------------------
tpm         = cfg_files;
tpm.tag     = 'tpm';
tpm.name    = 'Tissue probability map';
tpm.help    = {
               'Select the tissue probability image for this class. This should be map of eg grey matter, white matter or cerebro-spinal fluid probability. A nonlinear deformation field is estimated that best overlays the tissue probability maps on the individual subjects'' images. The default tissue probability maps for Unified Segmentation are modified versions of the ICBM Tissue Probabilistic Atlases.These tissue probability maps are kindly provided by the International Consortium for Brain Mapping, John C. Mazziotta and Arthur W. Toga. http://www.loni.ucla.edu/ICBM/ICBM_TissueProb.html. The original data are derived from 452 T1-weighted scans, which were aligned with an atlas space, corrected for scan inhomogeneities, and classified into grey matter, white matter and cerebrospinal fluid. These data were then affine registered to the MNI space and down-sampled to 2mm resolution. It is stored under SPM8_DIR/tpm/.'
                ''
                'The default TPM for New Segment has six tissues: grey matter, white matter, cerebro-spinal fluid, skull, scalp and air. It needs further refinement. The current versions were crudely generated (by JA) using data that was kindly provided by Cynthia Jongen of the Imaging Sciences Institute at Utrecht, NL. It is stored under SPM8_DIR/toolbox/Seg/.'
                ''
                'Also, an improved version of the TPM for New Segment is available at http://neuralengr.com/segment/. It has an extended field of view that covers the entire head (down to the jaw). It was generated by Chris Rorden at Georgia Institute of Technology, Atlanta, GA. A copy of this TPM is also available at SPM8_DIR/toolbox/mars/'
                ''
                'The model is refined further by allowing the tissue probability maps to be deformed according to a set of estimated parameters. This allows spatial normalisation and segmentation to be combined into the same model.'
               }';
tpm.filter = 'image';
tpm.ufilter = '.*';
tpm.num     = [1 1];
% ---------------------------------------------------------------------
% ngaus Num. Gaussians
% ---------------------------------------------------------------------
ngaus         = cfg_menu;
ngaus.tag     = 'ngaus';
ngaus.name    = 'Num. Gaussians';
ngaus.help    = {
                 'The number of Gaussians used to represent the intensity distribution for each tissue class can be greater than one. In other words, a tissue probability map may be shared by several clusters. The assumption of a single Gaussian distribution for each class does not hold for a number of reasons. In particular, a voxel may not be purely of one tissue type, and instead contain signal from a number of different tissues (partial volume effects). Some partial volume voxels could fall at the interface between different classes, or they may fall in the middle of structures such as the thalamus, which may be considered as being either grey or white matter. Various other image segmentation approaches use additional clusters to model such partial volume effects. These generally assume that a pure tissue class has a Gaussian intensity distribution, whereas intensity distributions for partial volume voxels are broader, falling between the intensities of the pure classes. Unlike these partial volume segmentation approaches, the model adopted here simply assumes that the intensity distribution of each class may not be Gaussian, and assigns belonging probabilities according to these non-Gaussian distributions. Typical numbers of Gaussians could be two for grey matter, two for white matter, two for CSF, three for bone, four for other soft tissues and two for air (background).'
                 ''
                 'The text above is excerpt from Ashburner and Friston 2005 Unified Segmentation paper. It is implemented in the New Segment toolbox, and also inherited here. However, this implementation is not corresponding to the proposed equation (Eq.27) in Ashburner and Friston 2005 paper. The correct equation corresponding to this implementation can be found at Eq.3c in Huang Y, Parra LC, Fully Automated Whole-Head Segmentation with Improved Smoothness and Continuity, with Theory Reviewed. PLoS ONE 2015, 10(5): e0125477. For implementation of Eq.27 in Ashburner and Friston 2005 paper, one is referred to SPM12 Segment toolbox.'
                 ''
                 'Note that if any of the Num. Gaussians is set to non-parametric, then a non-parametric approach will be used to model the tissue intensities. This may work for some images (eg CT), but not others - and it has not been optimised for multi-channel data. Note that it is likely to be especially problematic for images with poorly behaved intensity histograms due to aliasing effects that arise from having discrete values on the images.'
                 }';
ngaus.labels = {
                '1'
                '2'
                '3'
                '4'
                '5'
                '6'
                '7'
                '8'
                'Nonparametric'
                }';
ngaus.values = {
                1
                2
                3
                4
                5
                6
                7
                8
                Inf
                }';
ngaus.val    = {Inf};
% ---------------------------------------------------------------------
% native Native Tissue
% ---------------------------------------------------------------------
native         = cfg_menu;
native.tag     = 'native';
native.name    = 'Native Tissue';
native.help    = {'The native space option allows you to produce a tissue class image (c*) that is in alignment with the original/* (see Figure \ref{seg1})*/. It can also be used for ``importing'''' into a form that can be used with the DARTEL toolbox (rc*).'};
native.labels = {
                 'None'
                 'Native Space'
                 'DARTEL Imported'
                 'Native + DARTEL Imported'
                 }';
native.values = {
                 [0 0]
                 [1 0]
                 [0 1]
                 [1 1]
                 }';
native.val    = {[1 0]};
% ---------------------------------------------------------------------
% warped Warped Tissue
% ---------------------------------------------------------------------
warped         = cfg_menu;
warped.tag     = 'warped';
warped.name    = 'Warped Tissue';
warped.help    = {
                  'You can produce spatially normalised versions of the tissue class - both with (mwc*) and without (wc*) modulation (see below). These can be used for voxel-based morphometry. All you need to do is smooth them and do the stats.'
                  ''
                  '"Modulation" is to compensate for the effect of spatial normalisation.  When warping a series of images to match a template, it is inevitable that volumetric differences will be introduced into the warped images.  For example, if one subject''s temporal lobe has half the volume of that of the template, then its volume will be doubled during spatial normalisation. This will also result in a doubling of the voxels labelled grey matter.  In order to remove this confound, the spatially normalised grey matter (or other tissue class) is adjusted by multiplying by its relative volume before and after warping.  If warping results in a region doubling its volume, then the correction will halve the intensity of the tissue label. This whole procedure has the effect of preserving the total amount of grey matter signal in the normalised partitions.  Actually, in this version of SPM the warped data are not scaled by the Jacobian determinants when generating the "modulated" data.  Instead, the original voxels are projected into their new location in the warped images.  This exactly preserves the tissue count, but has the effect of introducing aliasing artifacts - especially if the original data are at a lower resolution than the warped images.  Smoothing should reduce this artifact though.'
                  ''
                  'Note also that the "unmodulated" data are generated slightly differently in this version of SPM. In this version, the projected data are corrected using a kind of smoothing procedure. This is not done exactly as it should be done (to save computational time), but it does a reasonable job. It also has the effect of extrapolating the warped tissue class images beyond the range of the original data.  This extrapolation is not perfect, as it is only an estimate, but it may still be a good thing to do.'
                  }';
warped.labels = {
                 'None'
                 'Modulated'
                 'Unmodulated'
                 'Modulated + Unmodulated'
                 }';
warped.values = {
                 [0 0]
                 [0 1]
                 [1 0]
                 [1 1]
                 }';
warped.val    = {[0 0]};
% ---------------------------------------------------------------------
% tissue Tissue
% ---------------------------------------------------------------------
tissue         = cfg_branch;
tissue.tag     = 'tissue';
tissue.name    = 'Tissue';
tissue.val     = {tpm ngaus native warped };
tissue.help    = {'A number of options are available for each of the tissues.  You may wish to save images of some tissues, but not others. If planning to use DARTEL, then make sure you generate ``imported'''' tissue class images of grey and white matter (and possibly others).  Different numbers of Gaussians may be needed to model the intensity distributions of the various tissues.'};
% ---------------------------------------------------------------------
% tissues Tissues
% ---------------------------------------------------------------------
tissues         = cfg_repeat;
tissues.tag     = 'tissues';
tissues.name    = 'Tissues';
tissues.values  = {tissue };
tissues.num     = [0 Inf];

tissues.val     = {tissue tissue tissue tissue tissue tissue };
tpm_nam = fullfile(fileparts(which(mfilename)),'eTPM.nii');
ngaus   = [2 2 2 3 4 2];
nval    = {[1 0],[1 0],[1 0],[1 0],[1 0],[1 0]};
for k=1:numel(ngaus),
    tissue.val{1}.val{1} = {[tpm_nam ',' num2str(k)]};
    tissue.val{2}.val    = {ngaus(k)};
    tissue.val{3}.val    = {nval{k}};
    tissues.val{k}       = tissue;
end

tissues.help = {'The data for each subject are classified into a number of different tissue types.  The tissue types are defined according to tissue probability maps, which define the prior probability of finding a tissue type at a particular location. Typically, the order of tissues is grey matter, white matter, CSF, bone, soft tissue and air/background (if using SPM8_DIR/toolbox/Seg/TPM.nii).'};

% % ---------------------------------------------------------------------
% % mrf MRF Parameter
% % ---------------------------------------------------------------------
% mrf         = cfg_entry;
% mrf.tag     = 'mrf';
% mrf.name    = 'MRF Parameter';
% mrf.help    = {'When tissue class images are written out, a few iterations of a simple Markov Random Field (MRF) cleanup procedure are run.  This parameter controls the strength of the MRF. Setting the value to zero will disable the cleanup.'};
% mrf.strtype = 'e';
% mrf.num     = [1 1];
% mrf.val     = {0}; % {2};
 % andy 2013-05-02

% % ---------------------------------------------------------------------
% % mrfIter Run MRF-New segment algorithm
% % ---------------------------------------------------------------------
% mrfIter         = cfg_entry;
% mrfIter.tag     = 'mrfIter';
% mrfIter.name    = 'MRF Iteration Tag';
% mrfIter.help    = {'This is a new algorithm for segmentation with MRF integrated based on the EM iterations using both TPM info and tissue correlation map (TCM). It uses the results from New Segment as initialization. Setting the value to one to use this algorithm, to zero to use default New Segment only.'};
% mrfIter.strtype = 'e';
% mrfIter.num     = [1 1];
% mrfIter.val     = {0}; % {1};
% % andy 2013-03-30
 % andy 2013-05-02

% ---------------------------------------------------------------------
% reg Warping Regularisation
% ---------------------------------------------------------------------
reg         = cfg_entry;
reg.tag     = 'reg';
reg.name    = 'Warping Regularisation';
reg.help    = {'The objective function for registering the tissue probability maps to the image to process, involves minimising the sum of two terms. One term gives a function of how probable the data is given the warping parameters. The other is a function of how probable the parameters are, and provides a penalty for unlikely deformations. Smoother deformations are deemed to be more probable. The amount of regularisation determines the tradeoff between the terms. Pick a value around one.  However, if your normalised images appear distorted, then it may be an idea to increase the amount of regularisation (by an order of magnitude). More regularisation gives smoother deformations, where the smoothness measure is determined by the bending energy of the deformations. '};
reg.strtype = 'e';
reg.num     = [1  1];
reg.val     = {4};
% ---------------------------------------------------------------------
% affreg Affine Regularisation
% ---------------------------------------------------------------------
affreg         = cfg_menu;
affreg.tag     = 'affreg';
affreg.name    = 'Affine Regularisation';
affreg.help    = {
                  'The procedure is a local optimisation, so it needs reasonable initial starting estimates. Images should be placed in approximate alignment using the Display function of SPM before beginning. A Mutual Information affine registration with the tissue probability maps (D''Agostino et al, 2004) is used to achieve approximate alignment. Note that this step does not include any model for intensity non-uniformity. This means that if the procedure is to be initialised with the affine registration, then the data should not be too corrupted with this artifact.If there is a lot of intensity non-uniformity, then manually position your image in order to achieve closer starting estimates, and turn off the affine registration.'
                  ''
                  'Affine registration into a standard space can be made more robust by regularisation (penalising excessive stretching or shrinking).  The best solutions can be obtained by knowing the approximate amount of stretching that is needed (e.g. ICBM templates are slightly bigger than typical brains, so greater zooms are likely to be needed). For example, if registering to an image in ICBM/MNI space, then choose this option.  If registering to a template that is close in size, then select the appropriate option for this.'
                  }';
affreg.labels = {
                 'No Affine Registration'
                 'ICBM space template - European brains'
                 'ICBM space template - East Asian brains'
                 'Average sized template'
                 'No regularisation'
                 }';
affreg.values = {
                 ''
                 'mni'
                 'eastern'
                 'subj'
                 'none'
                 }';
affreg.val    = {'mni'};
% ---------------------------------------------------------------------
% samp Sampling distance
% ---------------------------------------------------------------------
samp         = cfg_entry;
samp.tag     = 'samp';
samp.name    = 'Sampling distance';
samp.help    = {'This encodes the approximate distance between sampled points when estimating the model parameters. Smaller values use more of the data, but the procedure is slower and needs more memory. Determining the ``best'''' setting involves a compromise between speed and accuracy.'};
samp.strtype = 'e';
samp.num     = [1  1];
samp.val     = {3};
% ---------------------------------------------------------------------
% write Deformation Fields
% ---------------------------------------------------------------------
write         = cfg_menu;
write.tag     = 'write';
write.name    = 'Deformation Fields';
write.help    = {'Deformation fields can be saved to disk, and used by the Deformations Utility. For spatially normalising images to MNI space, you will need the forward deformation, whereas for spatially normalising (eg) GIFTI surface files, you''ll need the inverse. It is also possible to transform data in MNI space on to the individual subject, which also requires the inverse transform. Deformations are saved as .nii files, which contain three volumes to encode the x, y and z coordinates.'};
write.labels = {
                'None'
                'Inverse'
                'Forward'
                'Inverse + Forward'
                }';
write.values = {
                [0 0]
                [1 0]
                [0 1]
                [1 1]
                }';
write.val    = {[0 0]};
% ---------------------------------------------------------------------
% warp Warping
% ---------------------------------------------------------------------
warp         = cfg_branch;
warp.tag     = 'warp';
warp.name    = 'Warping';
warp.val     = {reg affreg samp write }; % andy 2013-05-02
warp.help    = {
'A number of warping options are provided, but the main one that you could consider changing is the one for specifying whether deformation fields or inverse deformation fields should be generated.'};

% ---------------------------------------------------------------------
% marsMode Switch to run different algorithms
% ---------------------------------------------------------------------
marsMode         = cfg_menu;
marsMode.tag     = 'marsMode';
marsMode.name    = 'MARS running mode';
marsMode.help    = {'A switch to run different algorithms. By default MARS will run. You can also just run the original New Segment provided by SPM8, or New Segment with MRF cleanup as a post-processing.'};
marsMode.labels  = {
                    'MARS'
                    'New Segment with MRF cleanup'
                    'New Segment'
                   }';
marsMode.values  = {
                    2
                    1
                    0
                   }';
marsMode.val     = {2};
 % andy 2015-07-16
 
% ---------------------------------------------------------------------
% tcm Tissue correlation map
% ---------------------------------------------------------------------
tcm         = cfg_files;
tcm.tag     = 'tcm';
tcm.name    = 'Tissue correlation map';
tcm.help    = {'Select the Tissue Correlation Map (TCM) for MARS iterations, or the MRF cleanup after New Segment. It can be a global TCM (K x K matrix), or local TCM (dim1 x dim2 x dim3 x K x K x 6 matrix), or regional TCM (K x K x 2 matrix). K is the number of tissue types in the segmentation results. If you need to run original New Segment which does not need any TCM, keep this option default, and the program will ignore it. If you need to run "New Segment with MRF cleanup", you can provide a global TCM as the MRF constraint for the MRF cleanup (local and regional TCMs are not supported), or you can just keep this option default and the program will use SPM8 default config to run "New Segment with MRF cleanup".'
               ''
               'Note that now MARS only supports 6-neighbourhood for the local TCM, i.e., right, left, anterior, posterior, superior and inferior neighbouring voxels. The dataset used to make the local TCM should be the same dataset that is used to make the TPM, otherwise the registration of the local TCM to the image that is to be segmented may go wrong. This is because the registration is estimated from the header information of TPM, and local TCM is assumed to have the same header info as the TPM and is registered by the same registration parameters estimated from TPM. For more details on how to make TPM/TCM from your own dataset, see Equation 4 in Huang Y, Parra LC, Fully Automated Whole-Head Segmentation with Improved Smoothness and Continuity, with Theory Reviewed. PLoS ONE 2015, 10(5): e0125477. Also note that the local TCM only works if it is oriented as RAS (right-anterior-superior). In other words, the heads that are used to generate local TCM and TPM should all be in RAS orientation (neurological convention). For global and regional TCM, there is no such restrictions.'};
tcm.filter = '.mat';
tcm.ufilter = '.*';
tcm.num     = [1 1];
% tcm.val     = {[]}; % or user load
 % andy 2015-07-16

% ---------------------------------------------------------------------
% beta The parameter to adjust the strength of local TCM
% ---------------------------------------------------------------------
beta         = cfg_entry;
beta.tag     = 'beta';
beta.name    = 'Beta for adjusting the strength of the local TCM';
beta.help    = {'Only applies to local TCM. Since it is intractable to obtain the analytical solutions of the energy terms from the TPM/TCM for the general N-site case, here in MARS we take a pragmatic approach, i.e., using a free parameter beta before interaction energy J_ij to balance the relative contribution of the neighbourhood prior versus the atlas. The default vaule of beta is 0.1. The higher beta is, the stronger effect the local TCM will have on the segmentation results, i.e., the segmentation will become more smooth. Avoid using higher beta to prevent over-regularisation.'};
beta.strtype = 'e';
beta.num     = [1 1];
beta.val     = {0.1};
 % andy 2015-07-16
 
% ---------------------------------------------------------------------
% marsConvergence User-defined convergence threshold for MARS iterations
% ---------------------------------------------------------------------
marsConvergence         = cfg_entry;
marsConvergence.tag     = 'marsConvergence';
marsConvergence.name    = 'User-defined convergence threshold for MARS iterations';
marsConvergence.help    = {'The convergence threshold for MARS iterations. It is defined as the maximal relative change of tissue volumes between consecutive iterations. Lower value means longer waiting time until convergence. Default value is set to 0.1, which only needs several minutes to convergence. Decrease the value to get better results upon convergence. A value of 1e-04 will take about 1 hour to convergence. Note here the estimated time does not include the time spent by New Segment.'};
marsConvergence.strtype = 'e';
marsConvergence.num     = [1 1];
marsConvergence.val     = {0.1};
 % andy 2015-07-16

% ---------------------------------------------------------------------
% marsOptions MARS related options
% ---------------------------------------------------------------------
marsOptions         = cfg_branch;
marsOptions.tag     = 'marsOptions';
marsOptions.name    = 'MARS related options';
marsOptions.val     = {marsMode tcm beta marsConvergence };
marsOptions.help    = {'These are the options specifically for MARS. You can specify whether to run MARS or just to run original New Segment in "MARS running mode" option.'};
 % andy 2015-07-16

marsOptions.val{2}.val{1} = {fullfile(fileparts(which(mfilename)),'rTCM_BW20_S1.mat')};
% ---------------------------------------------------------------------
% mars Morphologically and Anatomically accuRate Segmentation (MARS)
% ---------------------------------------------------------------------
mars         = cfg_exbranch;
mars.tag     = 'mars';
mars.name    = 'MARS';
mars.val     = {data tissues warp marsOptions};
mars.help    = {
                    'This toolbox is developed based on the SPM8 toolbox New Segment, which is an implementation of the Unified Segmentation algorithm (Ashburner J, Friston KJ. Unified segmentation. NeuroImage. 2005 Jul; 26(3):839–851). Here in this toolbox, the algorithm is further improved by adding the morphological constraints on neighbouring tissue voxels, which is encoded by Markov Random Field (MRF). This MRF information is added into the updating equation of the E-step of the algorithm, not just as a post-processing step. Compared to the New Segment, the resulting segmentation is shown to be more smooth and with less discontinuities, while at the same time the accuracy is not significantly changed. Therefore, it is named as Morphologically and Anatomically accuRate Segmentation (MARS). For more details on the algorithm, see Huang Y, Parra LC, Fully Automated Whole-Head Segmentation with Improved Smoothness and Continuity, with Theory Reviewed. PLoS ONE 2015, 10(5): e0125477.'
                    ''
                    'MARS is always initialised by the New Segment. In other words, this toolbox will run New Segment in the first step to get the optimal estimates on the bias field, the registration parameters between atlas and the image(s) to be segmented, and the initial estimates of segmentation posteriors and parameters of Gaussian mixture model. MARS will then continues to update the segmentation posteriors and Gaussian mixture model by using both atlas (tissue probability map, TPM) and tissue correlation map (TCM), while keep the bias field and registration parameters untouched. The toolbox is backwards compatible with the New Segment toolbox, as users can disable the MARS update so that only New Segment will be run.'
                    ''
                    'The TCM is a way of representing MRF constraint. It indicates the probability of co-occurrence of specific two tissues in two neighbouring voxels. For example, white matter cannot be adjacent to the air, and thus the probability of observing white matter as a neighbour of air voxel is 0, and this is encoded in the TCM. We developed three categories of TCM: the local, global and regional. The local TCM gives different co-occurrence probabilities depending on the location of the voxels, while the global version uses the same value for all the voxels in the image. The local TCM has nice theoretical properties, but not practical for implementation, as it needs too much memory to store and compute; while the global one is too rough to include any local information. The regional TCM is a trade-off between the two: it encodes locality information such as the interior of white matter and the boundary between grey and white matters, while at the same time uses the same co-occurrence value for one specific locality. Therefore, the regional TCM is highly recommended, and local TCM is not recommended unless you have a 64-bit computer with at least 50 GB memory.'
                    ''
                    'As this toolbox is directly inherited from the New Segment toolbox, it shares the same features, except the MARS update part. Some of the options in the toolbox do not yet work (which? Need to ask JA), and it has not yet been seamlessly integrated into the SPM8 software.'
                    }';
mars.prog = @spm_local_mars_run;
% mars.vout = @vout;
%----------------------------------------------------------------------

%======================================================================
function varargout = spm_local_mars_run(job)
varargout{1} = spm_mars_run(job);

% %======================================================================
% function dep = vout(job)
% % This depends on job contents, which may not be present when virtual
% % outputs are calculated.
% 
% cdep = cfg_dep;
% cdep(end).sname      = 'Seg Params';
% cdep(end).src_output = substruct('.','param','()',{':'});
% cdep(end).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
% 
% for i=1:numel(job.channel),
%     if job.channel(i).write(1),
%         cdep(end+1)          = cfg_dep;
%         cdep(end).sname      = sprintf('Bias Field (%d)',i);
%         cdep(end).src_output = substruct('.','channel','()',{i},'.','biasfield','()',{':'});
%         cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     end
%     if job.channel(i).write(2),
%         cdep(end+1)          = cfg_dep;
%         cdep(end).sname      = sprintf('Bias Corrected (%d)',i);
%         cdep(end).src_output = substruct('.','channel','()',{i},'.','biascorr','()',{':'});
%         cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     end
% end
% 
% for i=1:numel(job.tissue),
%     if job.tissue(i).native(1),
%         cdep(end+1)          = cfg_dep;
%         cdep(end).sname      = sprintf('c%d Images',i);
%         cdep(end).src_output = substruct('.','tiss','()',{i},'.','c','()',{':'});
%         cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     end
%     if job.tissue(i).native(2),
%         cdep(end+1)          = cfg_dep;
%         cdep(end).sname      = sprintf('rc%d Images',i);
%         cdep(end).src_output = substruct('.','tiss','()',{i},'.','rc','()',{':'});
%         cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     end
%     if job.tissue(i).warped(1),
%         cdep(end+1)          = cfg_dep;
%         cdep(end).sname      = sprintf('wc%d Images',i);
%         cdep(end).src_output = substruct('.','tiss','()',{i},'.','wc','()',{':'});
%         cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     end
%     if job.tissue(i).warped(2),
%         cdep(end+1)          = cfg_dep;
%         cdep(end).sname      = sprintf('mwc%d Images',i);
%         cdep(end).src_output = substruct('.','tiss','()',{i},'.','mwc','()',{':'});
%         cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     end
% end
% 
% if job.warp.write(1),
%     cdep(end+1)          = cfg_dep;
%     cdep(end).sname      = 'Inverse Deformations';
%     cdep(end).src_output = substruct('.','invdef','()',{':'});
%     cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
% end
% 
% if job.warp.write(2),
%     cdep(end+1)          = cfg_dep;
%     cdep(end).sname      = 'Forward Deformations';
%     cdep(end).src_output = substruct('.','fordef','()',{':'});
%     cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
% end
% 
% dep = cdep;
% 
