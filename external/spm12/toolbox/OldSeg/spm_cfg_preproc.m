function preproc = spm_cfg_preproc
% SPM Configuration file for toolbox 'Old Segment'
%______________________________________________________________________
% Copyright (C) 2005-2016 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_preproc.m 6894 2016-09-30 16:48:46Z spm $

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','OldSeg')); end

% ---------------------------------------------------------------------
% data Data
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'Data';
data.help    = {'Select scans for processing. This assumes that there is one scan for each subject. Note that multi-spectral (when there are two or more registered images of different contrasts) processing is not yet implemented for this method.'};
data.filter = 'image';
data.ufilter = '.*';
data.num     = [1 Inf];
% ---------------------------------------------------------------------
% GM Grey Matter
% ---------------------------------------------------------------------
GM         = cfg_menu;
GM.tag     = 'GM';
GM.name    = 'Grey Matter';
GM.help    = {'Options to produce grey matter images: c1*, wc1* and mwc1*.'};
GM.labels = {
             'None'
             'Native Space'
             'Unmodulated Normalised'
             'Modulated Normalised'
             'Native + Unmodulated Normalised'
             'Native + Modulated Normalised'
             'Native + Modulated + Unmodulated'
             'Modulated + Unmodulated Normalised'
}';
GM.values = {[0 0 0] [0 0 1] [0 1 0] [1 0 0] [0 1 1] [1 0 1] [1 1 1]...
             [1 1 0]};
GM.def    = @(val)spm_get_defaults('old.preproc.output.GM', val{:});
% ---------------------------------------------------------------------
% WM White Matter
% ---------------------------------------------------------------------
WM         = cfg_menu;
WM.tag     = 'WM';
WM.name    = 'White Matter';
WM.help    = {'Options to produce white matter images: c2*, wc2* and mwc2*.'};
WM.labels = {
             'None'
             'Native Space'
             'Unmodulated Normalised'
             'Modulated Normalised'
             'Native + Unmodulated Normalised'
             'Native + Modulated Normalised'
             'Native + Modulated + Unmodulated'
             'Modulated + Unmodulated Normalised'
}';
WM.values = {[0 0 0] [0 0 1] [0 1 0] [1 0 0] [0 1 1] [1 0 1] [1 1 1]...
             [1 1 0]};
WM.def    = @(val)spm_get_defaults('old.preproc.output.WM', val{:});
% ---------------------------------------------------------------------
% CSF Cerebro-Spinal Fluid
% ---------------------------------------------------------------------
CSF         = cfg_menu;
CSF.tag     = 'CSF';
CSF.name    = 'Cerebro-Spinal Fluid';
CSF.help    = {'Options to produce CSF images: c3*, wc3* and mwc3*.'};
CSF.labels = {
              'None'
              'Native Space'
              'Unmodulated Normalised'
              'Modulated Normalised'
              'Native + Unmodulated Normalised'
              'Native + Modulated Normalised'
              'Native + Modulated + Unmodulated'
              'Modulated + Unmodulated Normalised'
}';
CSF.values = {[0 0 0] [0 0 1] [0 1 0] [1 0 0] [0 1 1] [1 0 1] [1 1 1]...
             [1 1 0]};
CSF.def    = @(val)spm_get_defaults('old.preproc.output.CSF', val{:});
% ---------------------------------------------------------------------
% biascor Bias Corrected
% ---------------------------------------------------------------------
biascor         = cfg_menu;
biascor.tag     = 'biascor';
biascor.name    = 'Bias Corrected';
biascor.help    = {'This is the option to produce a bias corrected version of your image. MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.  The bias corrected version should have more uniform intensities within the different types of tissues.'};
biascor.labels = {
                  'Save Bias Corrected'
                  'Don''t Save Corrected'
}';
biascor.values = {1 0};
biascor.def    = @(val)spm_get_defaults('old.preproc.output.biascor', val{:});
% ---------------------------------------------------------------------
% cleanup Clean up any partitions
% ---------------------------------------------------------------------
cleanup         = cfg_menu;
cleanup.tag     = 'cleanup';
cleanup.name    = 'Clean up any partitions';
cleanup.help    = {
                   'This uses a crude routine for extracting the brain from segmented images.  It begins by taking the white matter, and eroding it acouple of times to get rid of any odd voxels.  The algorithm continues on to do conditional dilations for several iterations, where the condition is based upon gray or white matter being present. This identified region is then used to clean up the grey and white matter partitions, and has a slight influences on the CSF partition.'
                   ''
                   'If you find pieces of brain being chopped out in your data, then you may wish to disable or tone down the cleanup procedure.'
}';
cleanup.labels = {
                  'Dont do cleanup'
                  'Light Clean'
                  'Thorough Clean'
}';
cleanup.values = {0 1 2};
cleanup.def    = @(val)spm_get_defaults('old.preproc.output.cleanup', val{:});
% ---------------------------------------------------------------------
% output Output Files
% ---------------------------------------------------------------------
output         = cfg_branch;
output.tag     = 'output';
output.name    = 'Output Files';
output.val     = {GM WM CSF biascor cleanup };
output.help    = {
                  'This routine produces spatial normalisation parameters (*_seg_sn.mat files) by default. These can be used for writing spatially normalised versions of your data, via the "Normalise: Write" option. This mechanism may produce superior results than the "Normalise: Estimate" option (but probably not as good as those produced using DARTEL).'
                  ''
                  'In addition, it also produces files that can be used for doing inverse normalisation. If you have an image of regions defined in the standard space, then the inverse deformations can be used to warp these regions so that it approximately overlay your image. To use this facility, the bounding-box and voxel sizes should be set to non-finite values (e.g. [NaN NaN NaN] for the voxel sizes, and ones(2,3)*NaN for the bounding box. This would be done by the spatial normalisation module, which allows you to select a set of parameters that describe the nonlinear warps, and the images that they should be applied to.'
                  ''
                  'There are a number of options about what data you would like the routine to produce. The routine can be used for producing images of tissue classes, as well as bias corrected images. The native space option will produce a tissue class image (c*) that is in alignment with the original/* (see Figure \ref{seg1})*/.  You can also produce spatially normalised versions - both with (mwc*) and without (wc*) modulation/* (see Figure \ref{seg2})*/. The bounding box and voxel sizes of the spatially normalised versions are the same as that of the tissue probability maps with which they are registered. These can be used for doing voxel-based morphometry with (also see the ``Using DARTEL'' chapter of the manual). All you need to do is smooth them and do the stats (which means no more questions on the mailing list about how to do "optimized VBM").'
                  ''
                  'Modulation is to compensate for the effect of spatial normalisation.  When warping a series of images to match a template, it is inevitable that volumetric differences will be introduced into the warped images.  For example, if one subject''s temporal lobe has half the volume of that of the template, then its volume will be doubled during spatial normalisation. This will also result in a doubling of the voxels labelled grey matter.  In order to remove this confound, the spatially normalised grey matter (or other tissue class) is adjusted by multiplying by its relative volume before and after warping.  If warping results in a region doubling its volume, then the correction will halve the intensity of the tissue label. This whole procedure has the effect of preserving the total amount of grey matter signal in the normalised partitions.'
                  '/*\begin{figure} \begin{center} \includegraphics[width=140mm]{images/seg1} \end{center} \caption{Segmentation results. These are the results that can be obtained in the original space of the image (i.e. the results that are not spatially normalised). Top left: original image (X.img). Top right: bias corrected image (mX.img). Middle and bottom rows: segmented grey matter (c1X.img), white matter (c2X.img) and CSF (c3X.img). \label{seg1}} \end{figure} */'
                  '/*\begin{figure} \begin{center} \includegraphics[width=140mm]{images/seg2} \end{center} \caption{Segmentation results. These are the spatially normalised results that can be obtained (note that CSF data is not shown). Top row: The tissue probability maps used to guide the segmentation. Middle row: Spatially normalised tissue maps of grey and white matter (wc1X.img and wc2X.img). Bottom row: Modulated spatially normalised tissue maps of grey and white matter (mwc1X.img and mwc2X.img). \label{seg2}} \end{figure} */'
                  'A deformation field is a vector field, where three values are associated with each location in the field.  The field maps from co-ordinates in the normalised image back to co-ordinates in the original image.  The value of the field at co-ordinate [x y z] in the normalised space will be the co-ordinate [x'' y'' z''] in the original volume. The gradient of the deformation field at a co-ordinate is its Jacobian matrix, and it consists of a 3x3 matrix:'
                  ''
                  '%   /                      \'
                  '%   | dx''/dx  dx''/dy dx''/dz |'
                  '%   |                       |'
                  '%   | dy''/dx  dy''/dy dy''/dz |'
                  '%   |                       |'
                  '%   | dz''/dx  dz''/dy dz''/dz |'
                  '%   \                      /'
                  '/* \begin{eqnarray*}\begin{pmatrix}\frac{dx''}{dx} & \frac{dx''}{dy} & \frac{dx''}{dz}\cr\frac{dy''}{dx} & \frac{dy''}{dy} & \frac{dy''}{dz}\cr\frac{dz''}{dx} & \frac{dz''}{dy} & \frac{dz''}{dz}\cr\end{pmatrix}\end{eqnarray*}*/'
                  'The value of dx''/dy is a measure of how much x'' changes if y is changed by a tiny amount. The determinant of the Jacobian is the measure of relative volumes of warped and unwarped structures.  The modulation step simply involves multiplying by the relative volumes /*(see Figure \ref{seg2})*/.'
}';
% ---------------------------------------------------------------------
% tpm Tissue probability maps
% ---------------------------------------------------------------------
tpm         = cfg_files;
tpm.tag     = 'tpm';
tpm.name    = 'Tissue probability maps';
tpm.help    = {
               'Select the tissue probability images. These should be maps of grey matter, white matter and cerebro-spinal fluid probability. A nonlinear deformation field is estimated that best overlays the tissue probability maps on the individual subjects'' image. The default tissue probability maps are modified versions of the ICBM Tissue Probabilistic Atlases.These tissue probability maps are kindly provided by the International Consortium for Brain Mapping, John C. Mazziotta and Arthur W. Toga. http://www.loni.ucla.edu/ICBM/ICBM_TissueProb.html. The original data are derived from 452 T1-weighted scans, which were aligned with an atlas space, corrected for scan inhomogeneities, and classified into grey matter, white matter and cerebrospinal fluid. These data were then affine registered to the MNI space and downsampled to 2mm resolution.'
               ''
               'Rather than assuming stationary prior probabilities based upon mixing proportions, additional information is used, based on other subjects'' brain images.  Priors are usually generated by registering a large number of subjects together, assigning voxels to different tissue types and averaging tissue classes over subjects. Three tissue classes are used: grey matter, white matter and cerebro-spinal fluid. A fourth class is also used, which is simply one minus the sum of the first three. These maps give the prior probability of any voxel in a registered image being of any of the tissue classes - irrespective of its intensity.'
               ''
               'The model is refined further by allowing the tissue probability maps to be deformed according to a set of estimated parameters. This allows spatial normalisation and segmentation to be combined into the same model. This implementation uses a low-dimensional approach, which parameterises the deformations by a linear combination of about a thousand cosine transform bases. This is not an especially precise way of encoding deformations, but it can model the variability of overall brain shape. Evaluations by Hellier et al have shown that this simple model can achieve a registration accuracy comparable to other fully automated methods with many more parameters.'
}';
tpm.filter  = 'image';
tpm.dir     = fileparts(mfilename('fullpath'));
tpm.ufilter = '.*';
tpm.num     = [3 3];
tpm.def     = @(val)spm_get_defaults('old.preproc.tpm', val{:});
% ---------------------------------------------------------------------
% ngaus Gaussians per class
% ---------------------------------------------------------------------
ngaus         = cfg_entry;
ngaus.tag     = 'ngaus';
ngaus.name    = 'Gaussians per class';
ngaus.help    = {'The number of Gaussians used to represent the intensity distribution for each tissue class can be greater than one. In other words, a tissue probability map may be shared by several clusters. The assumption of a single Gaussian distribution for each class does not hold for a number of reasons. In particular, a voxel may not be purely of one tissue type, and instead contain signal from a number of different tissues (partial volume effects). Some partial volume voxels could fall at the interface between different classes, or they may fall in the middle of structures such as the thalamus, which may be considered as being either grey or white matter. Various other image segmentation approaches use additional clusters to model such partial volume effects. These generally assume that a pure tissue class has a Gaussian intensity distribution, whereas intensity distributions for partial volume voxels are broader, falling between the intensities of the pure classes. Unlike these partial volume segmentation approaches, the model adopted here simply assumes that the intensity distribution of each class may not be Gaussian, and assigns belonging probabilities according to these non-Gaussian distributions. Typical numbers of Gaussians could be two for grey matter, two for white matter, two for CSF, and four for everything else.'};
ngaus.strtype = 'n';
ngaus.num     = [4 1];
ngaus.def     = @(val)spm_get_defaults('old.preproc.ngaus', val{:});
% ---------------------------------------------------------------------
% regtype Affine Regularisation
% ---------------------------------------------------------------------
regtype         = cfg_menu;
regtype.tag     = 'regtype';
regtype.name    = 'Affine Regularisation';
regtype.help    = {
                   'The procedure is a local optimisation, so it needs reasonable initial starting estimates. Images should be placed in approximate alignment using the Display function of SPM before beginning. A Mutual Information affine registration with the tissue probability maps (D''Agostino et al, 2004) is used to achieve approximate alignment. Note that this step does not include any model for intensity non-uniformity. This means that if the procedure is to be initialised with the affine registration, then the data should not be too corrupted with this artifact.If there is a lot of intensity non-uniformity, then manually position your image in order to achieve closer starting estimates, and turn off the affine registration.'
                   ''
                   'Affine registration into a standard space can be made more robust by regularisation (penalising excessive stretching or shrinking).  The best solutions can be obtained by knowing the approximate amount of stretching that is needed (e.g. ICBM templates are slightly bigger than typical brains, so greater zooms are likely to be needed). For example, if registering to an image in ICBM/MNI space, then choose this option.  If registering to a template that is close in size, then select the appropriate option for this.'
}';
regtype.labels = {
                  'No Affine Registration'
                  'ICBM space template - European brains'
                  'ICBM space template - East Asian brains'
                  'Average sized template'
                  'No regularisation'
}';
regtype.values = {
                  ''
                  'mni'
                  'eastern'
                  'subj'
                  'none'
}';
regtype.def     = @(val)spm_get_defaults('old.preproc.regtype', val{:});
% ---------------------------------------------------------------------
% warpreg Warping Regularisation
% ---------------------------------------------------------------------
warpreg         = cfg_entry;
warpreg.tag     = 'warpreg';
warpreg.name    = 'Warping Regularisation';
warpreg.help    = {'The objective function for registering the tissue probability maps to the image to process, involves minimising the sum of two terms. One term gives a function of how probable the data is given the warping parameters. The other is a function of how probable the parameters are, and provides a penalty for unlikely deformations. Smoother deformations are deemed to be more probable. The amount of regularisation determines the tradeoff between the terms. Pick a value around one.  However, if your normalised images appear distorted, then it may be an idea to increase the amount of regularisation (by an order of magnitude). More regularisation gives smoother deformations, where the smoothness measure is determined by the bending energy of the deformations. '};
warpreg.strtype = 'e';
warpreg.num     = [1 1];
warpreg.def     = @(val)spm_get_defaults('old.preproc.warpreg', val{:});
% ---------------------------------------------------------------------
% warpco Warp Frequency Cutoff
% ---------------------------------------------------------------------
warpco         = cfg_entry;
warpco.tag     = 'warpco';
warpco.name    = 'Warp Frequency Cutoff';
warpco.help    = {'Cutoff of DCT bases.  Only DCT bases of periods longer than the cutoff are used to describe the warps. The number actually used will depend on the cutoff and the field of view of your image. A smaller cutoff frequency will allow more detailed deformations to be modelled, but unfortunately comes at a cost of greatly increasing the amount of memory needed, and the time taken.'};
warpco.strtype = 'e';
warpco.num     = [1 1];
warpco.def     = @(val)spm_get_defaults('old.preproc.warpco', val{:});
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
biasreg.values = {0 1e-5 .0001 .001 .01 .1 1 10};
biasreg.def     = @(val)spm_get_defaults('old.preproc.biasreg', val{:});
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
biasfwhm.values = {30 40 50 60 70 80 90 100 110 120 130 140 150 Inf};
biasfwhm.def     = @(val)spm_get_defaults('old.preproc.biasfwhm', val{:});
% ---------------------------------------------------------------------
% samp Sampling distance
% ---------------------------------------------------------------------
samp         = cfg_entry;
samp.tag     = 'samp';
samp.name    = 'Sampling distance';
samp.help    = {'The approximate distance between sampled points when estimating the model parameters. Smaller values use more of the data, but the procedure is slower.'};
samp.strtype = 'e';
samp.num     = [1 1];
samp.def     = @(val)spm_get_defaults('old.preproc.samp', val{:});
% ---------------------------------------------------------------------
% msk Masking image
% ---------------------------------------------------------------------
msk         = cfg_files;
msk.tag     = 'msk';
msk.name    = 'Masking image';
msk.val{1} = {''};
msk.help    = {'The segmentation can be masked by an image that conforms to the same space as the images to be segmented.  If an image is selected, then it must match the image(s) voxel-for voxel, and have the same voxel-to-world mapping.  Regions containing a value of zero in this image do not contribute when estimating the various parameters. '};
msk.filter = 'image';
msk.ufilter = '.*';
msk.num     = [0 1];
% ---------------------------------------------------------------------
% opts Custom
% ---------------------------------------------------------------------
opts         = cfg_branch;
opts.tag     = 'opts';
opts.name    = 'Custom';
opts.val     = {tpm ngaus regtype warpreg warpco biasreg biasfwhm samp msk };
opts.help    = {'Various options can be adjusted in order to improve the performance of the algorithm with your data.  Knowing what works best should be a matter of empirical exploration.  For example, if your data has very little intensity non-uniformity artifact, then the bias regularisation should be increased.  This effectively tells the algorithm that there is very little bias in your data, so it does not try to model it.'};
% ---------------------------------------------------------------------
% oldseg Old Segment
% ---------------------------------------------------------------------
preproc         = cfg_exbranch;
preproc.tag     = 'oldseg';
preproc.name    = 'Old Segment';
preproc.val     = {data output opts};
preproc.help    = {
                   'Segment, bias correct and spatially normalise - all in the same model/* \cite{ashburner05}*/. This function can be used for bias correcting, spatially normalising or segmenting your data. Note that this module needs the images to be roughly aligned with the tissue probability maps before you begin.  If strange results are obtained, then this is usually because the images were poorly aligned beforehand.  The Display option can be used to manually reposition the images so that the AC is close to coordinate 0,0,0 (within a couple of cm) and the orientation is within a few degrees of the tissue probability map data.'
                   ''
                   'Many investigators use tools within older versions of SPM for a technique that has become known as "optimised" voxel-based morphometry (VBM). VBM performs region-wise volumetric comparisons among populations of subjects. It requires the images to be spatially normalised, segmented into different tissue classes, and smoothed, prior to performing statistical tests/* \cite{wright_vbm,am_vbmreview,ashburner00b,john_should}*/. The "optimised" pre-processing strategy involved spatially normalising subjects'' brain images to a standard space, by matching grey matter in these images, to a grey matter reference.  The historical motivation behind this approach was to reduce the confounding effects of non-brain (e.g. scalp) structural variability on the registration. Tissue classification in older versions of SPM required the images to be registered with tissue probability maps. After registration, these maps represented the prior probability of different tissue classes being found at each location in an image.  Bayes rule can then be used to combine these priors with tissue type probabilities derived from voxel intensities, to provide the posterior probability.'
                   ''
                   'This procedure was inherently circular, because the registration required an initial tissue classification, and the tissue classification requires an initial registration.  This circularity is resolved here by combining both components into a single generative model. This model also includes parameters that account for image intensity non-uniformity. Estimating the model parameters (for a maximum a posteriori solution) involves alternating among classification, bias correction and registration steps. This approach provides better results than simple serial applications of each component.'
                   ''
                   'Note that multi-spectral segmentation (e.g. from a registered T1 and T2 image) is not yet implemented, but is planned for a future SPM version.'
}';
preproc.prog = @spm_run_preproc;
preproc.vout = @vout;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function dep = vout(job)
opts  = job.output;

cdep(1)            = cfg_dep;
cdep(1).sname      = 'Norm Params Subj->MNI';
cdep(1).src_output = substruct('()',{1}, '.','snfile','()',{':'});
cdep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
cdep(2)            = cfg_dep;
cdep(2).sname      = 'Norm Params MNI->Subj';
cdep(2).src_output = substruct('()',{1}, '.','isnfile','()',{':'});
cdep(2).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
if opts.biascor,
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Bias Corr Images';
    cdep(end).src_output = substruct('()',{1}, '.','biascorr','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end;
sonames = {'GM','WM','CSF'};
for k1=1:3,
    if ~strcmp(opts.(sonames{k1}),'<UNDEFINED>')
        if opts.(sonames{k1})(3),
            cdep(end+1)          = cfg_dep;
            cdep(end).sname      = sprintf('c%d Images',k1);
            cdep(end).src_output = substruct('()',{1}, '.',sprintf('c%d',k1),'()',{':'});
            cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        end;
        if opts.(sonames{k1})(2),
            cdep(end+1)          = cfg_dep;
            cdep(end).sname      = sprintf('wc%d Images',k1);
            cdep(end).src_output = substruct('()',{1}, '.',sprintf('wc%d',k1),'()',{':'});
            cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        end;
        if opts.(sonames{k1})(1),
            cdep(end+1)          = cfg_dep;
            cdep(end).sname      = sprintf('mwc%d Images',k1);
            cdep(end).src_output = substruct('()',{1}, '.',sprintf('mwc%d',k1),'()',{':'});
            cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        end;
    end
end;
dep = cdep;
