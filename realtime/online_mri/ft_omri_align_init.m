function model = ft_omri_align_init(Vr,flags)
% function model = ft_omri_align_init(Vr,flags)
%
% Ripped out of SPM 8 and modified (2010, S.Klanke)
%
% Estimation of within modality rigid body movement parameters
% FORMAT P = spm_realign(P,flags)
%
% P     - matrix of filenames {one string per row}
%         All operations are performed relative to the first image.
%         ie. Coregistration is to the first image, and resampling
%         of images is into the space of the first image.
%         For multiple sessions, P should be a cell array, where each
%         cell should be a matrix of filenames.
%
% flags - a structure containing various options.  The fields are:
%         quality - Quality versus speed trade-off.  Highest quality
%                   (1) gives most precise results, whereas lower
%                   qualities gives faster realignment.
%                   The idea is that some voxels contribute little to
%                   the estimation of the realignment parameters.
%                   This parameter is involved in selecting the number
%                   of voxels that are used.
%
%         fwhm    - The FWHM of the Gaussian smoothing kernel (mm)
%                   applied to the images before estimating the
%                   realignment parameters.
%
%         sep     - the default separation (mm) to sample the images.
%
%         PW      - a filename of a weighting image (reciprocal of
%                   standard deviation).  If field does not exist, then
%                   no weighting is done.
%
%         interp  - B-spline degree used for interpolation
%
%__________________________________________________________________________
%
% Inputs
% A series of *.img conforming to SPM data format (see 'Data Format').
%
% Outputs
% If no output argument, then an updated voxel to world matrix is written
% to the headers of the images (a .mat file is created for 4D images).
% The details of the transformation are displayed in the
% results window as plots of translation and rotation.
% A set of realignment parameters are saved for each session, named:
% rp_*.txt.
%__________________________________________________________________________
%
% The voxel to world mappings.
%
% These are simply 4x4 affine transformation matrices represented in the
% NIFTI headers (see http://nifti.nimh.nih.gov/nifti-1 ).
% These are normally modified by the `realignment' and `coregistration'
% modules.  What these matrixes represent is a mapping from
% the voxel coordinates (x0,y0,z0) (where the first voxel is at coordinate
% (1,1,1)), to coordinates in millimeters (x1,y1,z1).
%  
% x1 = M(1,1)*x0 + M(1,2)*y0 + M(1,3)*z0 + M(1,4)
% y1 = M(2,1)*x0 + M(2,2)*y0 + M(2,3)*z0 + M(2,4)
% z1 = M(3,1)*x0 + M(3,2)*y0 + M(3,3)*z0 + M(3,4)
%
% Assuming that image1 has a transformation matrix M1, and image2 has a
% transformation matrix M2, the mapping from image1 to image2 is: M2\M1
% (ie. from the coordinate system of image1 into millimeters, followed
% by a mapping from millimeters into the space of image2).
%
% These matrices allow several realignment or coregistration steps to be
% combined into a single operation (without the necessity of resampling the
% images several times).  The `.mat' files are also used by the spatial
% normalisation module.
%__________________________________________________________________________
% Ref:
% Friston KJ, Ashburner J, Frith CD, Poline J-B, Heather JD & Frackowiak
% RSJ (1995) Spatial registration and normalization of images Hum. Brain
% Map. 2:165-189
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$


if nargin==0, return; end;

model = struct('quality',1,'fwhm',5,'sep',4,'interp',2,'wrap',[0 0 0],'rtm',0,'PW','','lkp',1:6, 'mat', eye(4), 'time', 2);
if nargin > 1,
    fnms = fieldnames(model);
    for i=1:length(fnms),
        if isfield(flags,fnms{i}),
            model.(fnms{i}) = flags.(fnms{i});
        end;
    end;
end;


% assume reference volume is not scaled or anything
skip = sqrt(sum(model.mat(1:3,1:3).^2)).^(-1)*model.sep;
d    = size(Vr);
rand('state',0); % want the results to be consistant.
if d(3) < 3,
    model.lkp = [1 2 6];
    [x1,x2,x3] = ndgrid(1:skip(1):d(1)-.5, 1:skip(2):d(2)-.5, 1:skip(3):d(3));
    x1   = x1 + rand(size(x1))*0.5;
    x2   = x2 + rand(size(x2))*0.5;
else
    [x1,x2,x3]=ndgrid(1:skip(1):d(1)-.5, 1:skip(2):d(2)-.5, 1:skip(3):d(3)-.5);
    x1   = x1 + rand(size(x1))*0.5;
    x2   = x2 + rand(size(x2))*0.5;
    x3   = x3 + rand(size(x3))*0.5; 
end;

x1   = x1(:);
x2   = x2(:);
x3   = x3(:);

% Possibly mask an area of the sample volume.
%-----------------------------------------------------------------------
if ~isempty(model.PW),
	[y1,y2,y3]=coords([0 0 0  0 0 0],model.mat, model.PW.mat ,x1,x2,x3);
    wt  = spm_sample_vol(model.PW,y1,y2,y3,1);
    msk = find(wt>0.01);
    x1  = x1(msk);
    x2  = x2(msk);
    x3  = x3(msk);
    wt  = wt(msk);
else
    wt = [];
end;

% Compute rate of change of chi2 w.r.t changes in parameters (matrix A)
%-----------------------------------------------------------------------
V   = smooth_vol(Vr, model.interp, model.wrap, model.fwhm, model.mat);
deg = [model.interp*[1 1 1]' model.wrap(:)];

[G,dG1,dG2,dG3] = spm_bsplins(V,x1,x2,x3,deg);
clear V
A0 = make_A(model.mat,x1,x2,x3,dG1,dG2,dG3,wt,model.lkp);

b  = G;
if ~isempty(wt), b = b.*wt; end;

    % Remove voxels that contribute very little to the final estimate.
    % Simulated annealing or something similar could be used to
    % eliminate a better choice of voxels - but this way will do for
    % now. It basically involves removing the voxels that contribute
    % least to the determinant of the inverse covariance matrix.

%    spm_chi2_plot('Init','Eliminating Unimportant Voxels',...
              %'Relative quality','Iteration');
Alpha = [A0 b];
Alpha = Alpha'*Alpha;
det0  = det(Alpha);
det1  = det0;
%    spm_chi2_plot('Set',det1/det0);
while det1/det0 > model.quality,
	dets  = zeros(size(A0,1),1);
	for i=1:size(A0,1),
	   tmp     = [A0(i,:) b(i)];
		dets(i) = det(Alpha - tmp'*tmp);
	end;
	clear tmp
	[junk,msk] = sort(det1-dets);
	msk        = msk(1:round(length(dets)/10));
	A0(msk,:) = [];   b(msk,:) = [];   G(msk,:) = [];
	x1(msk,:) = [];  x2(msk,:) = [];  x3(msk,:) = [];
	dG1(msk,:) = []; dG2(msk,:) = []; dG3(msk,:) = [];
	if ~isempty(wt),  wt(msk,:) = []; end;
	Alpha = [A0 b];
	Alpha = Alpha'*Alpha;
	det1  = det(Alpha);
%	spm_chi2_plot('Set',single(det1/det0));
end;

model.deg = deg;
model.x1 = x1;
model.x2 = x2;
model.x3 = x3;
model.wt = wt;
model.A0 = A0;
model.b  = b;
model.mask = true(size(Vr));
