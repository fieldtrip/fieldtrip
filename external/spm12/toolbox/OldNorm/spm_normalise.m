function params = spm_normalise(VG,VF,matname,VWG,VWF,flags)
% Spatial (stereotactic) normalization
%
% FORMAT params = spm_normalise(VG,VF,matname,VWG,VWF,flags)
% VG        - template handle(s)
% VF        - handle of image to estimate params from
% matname   - name of file to store deformation definitions
% VWG       - template weighting image
% VWF       - source weighting image
% flags     - flags.  If any field is not passed, then defaults are assumed.
%             smosrc - smoothing of source image (FWHM of Gaussian in mm).
%                      Defaults to 8.
%             smoref - smoothing of template image (defaults to 0).
%             regtype - regularisation type for affine registration
%                       See spm_affreg.m (default = 'mni').
%             cutoff  - Cutoff of the DCT bases.  Lower values mean more
%                       basis functions are used (default = 30mm).
%             nits    - number of nonlinear iterations (default=16).
%             reg     - amount of regularisation (default=0.1)
% _________________________________________________________________________
% 
% This module spatially (stereotactically) normalizes MRI, PET or SPECT
% images into a standard space defined by some ideal model or template
% image[s].  The template images supplied with SPM conform to the space
% defined by the ICBM, NIH P-20 project, and approximate that of the
% the space described in the atlas of Talairach and Tournoux (1988).
% The transformation can also be applied to any other image that has
% been coregistered with these scans.
%
% 
% Mechanism
% Generally, the algorithms work by minimising the sum of squares
% difference between the image which is to be normalised, and a linear
% combination of one or more template images.  For the least squares
% registration to produce an unbiased estimate of the spatial
% transformation, the image contrast in the templates (or linear
% combination of templates) should be similar to that of the image from
% which the spatial normalization is derived.  The registration simply
% searches for an optimum solution.  If the starting estimates are not
% good, then the optimum it finds may not find the global optimum.
% 
% The first step of the normalization is to determine the optimum
% 12-parameter affine transformation.  Initially, the registration is
% performed by matching the whole of the head (including the scalp) to
% the template.  Following this, the registration proceeded by only
% matching the brains together, by appropriate weighting of the template
% voxels.  This is a completely automated procedure (that does not
% require ``scalp editing'') that discounts the confounding effects of
% skull and scalp differences.   A Bayesian framework is used, such that
% the registration searches for the solution that maximizes the a
% posteriori probability of it being correct.  i.e., it maximizes the
% product of the likelihood function (derived from the residual squared
% difference) and the prior function (which is based on the probability
% of obtaining a particular set of zooms and shears).
% 
% The affine registration is followed by estimating nonlinear deformations,
% whereby the deformations are defined by a linear combination of three
% dimensional discrete cosine transform (DCT) basis functions.
% The parameters represent coefficients of the deformations in
% three orthogonal directions.  The matching involved simultaneously
% minimizing the bending energies of the deformation fields and the
% residual squared difference between the images and template(s).
% 
% An option is provided for allowing weighting images (consisting of pixel
% values between the range of zero to one) to be used for registering
% abnormal or lesioned brains.  These images should match the dimensions
% of the image from which the parameters are estimated, and should contain
% zeros corresponding to regions of abnormal tissue.
% 
% 
% Uses
% Primarily for stereotactic normalization to facilitate inter-subject
% averaging and precise characterization of functional anatomy.  It is
% not necessary to spatially normalise the data (this is only a
% pre-requisite  for  intersubject averaging or reporting in the
% Talairach space).
% 
% Inputs
% The first input is the image which is to be normalised. This image
% should be of the same modality (and MRI sequence etc) as the template
% which is specified. The same spatial transformation can then be
% applied to any other images of the same subject.  These files should
% conform to the SPM data format (See 'Data Format'). Many subjects can
% be entered at once, and there is no restriction on image dimensions
% or voxel size.
% 
% Providing that the images have a correct ".mat" file associated with
% them, which describes the spatial relationship between them, it is
% possible to spatially normalise the images without having first
% resliced them all into the same space. The ".mat" files are generated
% by "spm_realign" or "spm_coregister".
% 
% Default values of parameters pertaining to the extent and sampling of
% the standard space can be changed, including the model or template
% image[s].
% 
% 
% Outputs
% All normalized *.img scans are written to the same subdirectory as
% the original *.img, prefixed with a 'n' (i.e. n*.img).  The details
% of the transformations are displayed in the results window, and the
% parameters are saved in the "*_sn.mat" file.
% 
%__________________________________________________________________________
%
% References:
% K.J. Friston, J. Ashburner, C.D. Frith, J.-B. Poline, J.D. Heather,
% and R.S.J. Frackowiak
% Spatial Registration and Normalization of Images.
% Human Brain Mapping 2:165-189, 1995.
% 
% J. Ashburner, P. Neelin, D.L. Collins, A.C. Evans and K.J. Friston
% Incorporating Prior Knowledge into Image Registration.
% NeuroImage 6:344-352, 1997.
%
% J. Ashburner and K.J. Friston
% Nonlinear spatial normalization using basis functions.
% Human Brain Mapping, 7(4):254-266, 1999.
%__________________________________________________________________________
% Copyright (C) 2002-2013 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_normalise.m 5745 2013-11-14 12:46:49Z guillaume $


if nargin<2, error('Incorrect usage.'); end;
if ischar(VF), VF = spm_vol(VF); end;
if ischar(VG), VG = spm_vol(VG); end;
if nargin<3,
    if nargout==0,
        [pth,nm]  = spm_fileparts(deblank(VF(1).fname));
        matname   = fullfile(pth,[nm '_sn.mat']);
    else
        matname   = '';
    end;
end;
if nargin<4, VWG  = ''; end;
if nargin<5, VWF  = ''; end;
if ischar(VWG), VWG=spm_vol(VWG); end;
if ischar(VWF), VWF=spm_vol(VWF); end;                                                                     


def_flags          = spm_get_defaults('old.normalise.estimate');
def_flags.graphics = 1;
if nargin < 6,
    flags = def_flags;
else
    fnms  = fieldnames(def_flags);
    for i=1:length(fnms),
        if ~isfield(flags,fnms{i}),
            flags.(fnms{i}) = def_flags.(fnms{i});
        end;
    end;
end;

fprintf('Smoothing by %g & %gmm..\n', flags.smoref, flags.smosrc);
VF1 = spm_smoothto8bit(VF,flags.smosrc);

% Rescale images so that globals are better conditioned
VF1.pinfo(1:2,:) = VF1.pinfo(1:2,:)/spm_global(VF1);
for i=1:numel(VG),
    VG1(i) = spm_smoothto8bit(VG(i),flags.smoref);
    VG1(i).pinfo(1:2,:) = VG1(i).pinfo(1:2,:)/spm_global(VG(i));
end;

% Affine Normalisation
%--------------------------------------------------------------------------
fprintf('Coarse Affine Registration..\n');
aflags    = struct('sep',max(flags.smoref,flags.smosrc), 'regtype',flags.regtype,...
    'WG',[],'WF',[],'globnorm',0);
aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));

M         = eye(4); %spm_matrix(prms');
spm_plot_convergence('Init','Affine Registration','Mean squared difference','Iteration');
[M,scal]  = spm_affreg(VG1, VF1, aflags, M);
 
fprintf('Fine Affine Registration..\n');
aflags.WG  = VWG;
aflags.WF  = VWF;
aflags.sep = aflags.sep/2;
[M,scal]   = spm_affreg(VG1, VF1, aflags, M,scal);
Affine     = inv(VG(1).mat\M*VF1(1).mat);
spm_plot_convergence('Clear');

% Basis function Normalisation
%--------------------------------------------------------------------------
fov = VF1(1).dim(1:3).*sqrt(sum(VF1(1).mat(1:3,1:3).^2));
if any(fov<15*flags.smosrc/2 & VF1(1).dim(1:3)<15),
    fprintf('Field of view too small for nonlinear registration\n');
    Tr = [];
elseif isfinite(flags.cutoff) && flags.nits && ~isinf(flags.reg),
        fprintf('3D CT Norm...\n');
    Tr = snbasis(VG1,VF1,VWG,VWF,Affine,...
        max(flags.smoref,flags.smosrc),flags.cutoff,flags.nits,flags.reg);
else
    Tr = [];
end;
clear VF1 VG1

flags.version = '$Rev: 5745 $';
flags.date    = date;

params = struct('Affine',Affine, 'Tr',Tr, 'VF',VF, 'VG',VG, 'flags',flags);

if flags.graphics, spm_normalise_disp(params,VF); end;

% Remove dat fields before saving
%--------------------------------------------------------------------------
if isfield(VF,'dat'), VF = rmfield(VF,'dat'); end;
if isfield(VG,'dat'), VG = rmfield(VG,'dat'); end;
if ~isempty(matname),
    fprintf('Saving Parameters..\n');
    save(matname,'Affine','Tr','VF','VG','flags', spm_get_defaults('mat.format'));
end;
return;


%==========================================================================
% function Tr = snbasis(VG,VF,VWG,VWF,Affine,fwhm,cutoff,nits,reg)
%==========================================================================
function Tr = snbasis(VG,VF,VWG,VWF,Affine,fwhm,cutoff,nits,reg)
% 3D Basis Function Normalization
% FORMAT Tr = snbasis(VG,VF,VWG,VWF,Affine,fwhm,cutoff,nits,reg)
% VG        - Template volumes (see spm_vol).
% VF        - Volume to normalize.
% VWG       - weighting Volume - for template.
% VWF       - weighting Volume - for object.
% Affine    - A 4x4 transformation (in voxel space).
% fwhm      - smoothness of images.
% cutoff    - frequency cutoff of basis functions.
% nits      - number of iterations.
% reg       - regularisation.
% Tr - Discrete cosine transform of the warps in X, Y & Z.
%
% snbasis performs a spatial normalization based upon a 3D
% discrete cosine transform.
%__________________________________________________________________________

fwhm    = [fwhm 30];

% Number of basis functions for x, y & z
%--------------------------------------------------------------------------
tmp  = sqrt(sum(VG(1).mat(1:3,1:3).^2));
k    = max(round((VG(1).dim(1:3).*tmp)/cutoff),[1 1 1]);

% Scaling is to improve stability.
%--------------------------------------------------------------------------
stabilise = 8;
basX = spm_dctmtx(VG(1).dim(1),k(1))*stabilise;
basY = spm_dctmtx(VG(1).dim(2),k(2))*stabilise;
basZ = spm_dctmtx(VG(1).dim(3),k(3))*stabilise;

dbasX = spm_dctmtx(VG(1).dim(1),k(1),'diff')*stabilise;
dbasY = spm_dctmtx(VG(1).dim(2),k(2),'diff')*stabilise;
dbasZ = spm_dctmtx(VG(1).dim(3),k(3),'diff')*stabilise;

vx1 = sqrt(sum(VG(1).mat(1:3,1:3).^2));
vx2 = vx1;
kx = (pi*((1:k(1))'-1)/VG(1).dim(1)/vx1(1)).^2; ox=ones(k(1),1);
ky = (pi*((1:k(2))'-1)/VG(1).dim(2)/vx1(2)).^2; oy=ones(k(2),1);
kz = (pi*((1:k(3))'-1)/VG(1).dim(3)/vx1(3)).^2; oz=ones(k(3),1);

if 1,
    % BENDING ENERGY REGULARIZATION
    % Estimate a suitable sparse diagonal inverse covariance matrix for
    % the parameters (IC0).
    %----------------------------------------------------------------------
    IC0 = (1*kron(kz.^2,kron(ky.^0,kx.^0)) +...
        1*kron(kz.^0,kron(ky.^2,kx.^0)) +...
        1*kron(kz.^0,kron(ky.^0,kx.^2)) +...
        2*kron(kz.^1,kron(ky.^1,kx.^0)) +...
        2*kron(kz.^1,kron(ky.^0,kx.^1)) +...
        2*kron(kz.^0,kron(ky.^1,kx.^1)) );
    IC0 = reg*IC0*stabilise^6;
    IC0 = [IC0*vx2(1)^4 ; IC0*vx2(2)^4 ; IC0*vx2(3)^4 ; zeros(prod(size(VG))*4,1)];
    IC0 = sparse(1:length(IC0),1:length(IC0),IC0,length(IC0),length(IC0));
else
    % MEMBRANE ENERGY (LAPLACIAN) REGULARIZATION
    %----------------------------------------------------------------------
    IC0 = kron(kron(oz,oy),kx) + kron(kron(oz,ky),ox) + kron(kron(kz,oy),ox);
    IC0 = reg*IC0*stabilise^6;
    IC0 = [IC0*vx2(1)^2 ; IC0*vx2(2)^2 ; IC0*vx2(3)^2 ; zeros(prod(size(VG))*4,1)];
    IC0 = sparse(1:length(IC0),1:length(IC0),IC0,length(IC0),length(IC0));
end;

% Generate starting estimates.
%--------------------------------------------------------------------------
s1 = 3*prod(k);
s2 = s1 + numel(VG)*4;
T  = zeros(s2,1);
T(s1+(1:4:numel(VG)*4)) = 1;

pVar = Inf;
for iter=1:nits,
    fprintf(' iteration %2d: ', iter);
    [Alpha,Beta,Var,fw] = spm_brainwarp(VG,VF,Affine,basX,basY,basZ,dbasX,dbasY,dbasZ,T,fwhm,VWG, VWF);
    if Var>pVar, scal = pVar/Var ; Var = pVar; else scal = 1; end;
    pVar = Var;
    T = (Alpha + IC0*scal)\(Alpha*T + Beta);
    fwhm(2) = min([fw fwhm(2)]);
    fprintf(' FWHM = %6.4g Var = %g\n', fw,Var);
end;

% Values of the 3D-DCT - for some bizarre reason, this needs to be done
% as two separate statements in Matlab 6.5...
%--------------------------------------------------------------------------
Tr = reshape(T(1:s1),[k 3]);
drawnow;
Tr = Tr*stabilise.^3;
return;
