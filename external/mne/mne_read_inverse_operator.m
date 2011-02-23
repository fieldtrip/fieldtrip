function [inv] = mne_read_inverse_operator(fname)
%
% [inv] = mne_read_inverse_operator(fname)
%
% Reads the inverse operator decomposition from a fif file
%
% fname        - The name of the file
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.8  2008/05/26 10:49:26  msh
%   Update to incorporate the already weighted lead field basis
%
%   Revision 1.7  2007/01/29 21:21:22  msh
%   Added reading of the additional source prior covariances.
%
%   Revision 1.6  2006/09/28 01:06:07  msh
%   Added nchan field to the inverse operator structure.
%
%   Revision 1.5  2006/05/05 03:50:40  msh
%   Added routines to compute L2-norm inverse solutions.
%   Added mne_write_inverse_sol_stc to write them in stc files
%   Several bug fixes in other files
%
%   Revision 1.4  2006/05/03 19:09:03  msh
%   Fixed even more compatibility issues.
%
%   Revision 1.3  2006/04/23 15:29:41  msh
%   Added MGH to the copyright
%
%   Revision 1.2  2006/04/21 21:33:07  msh
%   Fixed error in combination of the forward solution.
%   Fixed help text in mne_read_inverse_operator
%
%   Revision 1.1  2006/04/20 21:49:38  msh
%   Added mne_read_inverse_operator
%   Changed some of the routines accordingly for more flexibility.
%
%

me='MNE:mne_read_inverse_operator';

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end

if nargin ~= 1
    error(me,'Incorrect number of arguments');
end
%
%   Open the file, create directory
%
fprintf(1,'Reading inverse operator decomposition from %s...\n',fname);
[ fid, tree ] = fiff_open(fname);
%
%   Find all inverse operators
%
invs = fiff_dir_tree_find(tree,FIFF.FIFFB_MNE_INVERSE_SOLUTION);
if isempty(invs)
    fclose(fid);
    error(me,'No inverse solutions in %s',fname);
end
invs = invs(1);
%
%   Parent MRI data
%
parent_mri = fiff_dir_tree_find(tree,FIFF.FIFFB_MNE_PARENT_MRI_FILE);
if isempty(parent_mri)
    fclose(fid);
    error(me,'No parent MRI information in %s',fname);
end
fprintf(1,'\tReading inverse operator info...');
%
%   Methods and source orientations
%
tag = find_tag(invs,FIFF.FIFF_MNE_INCLUDED_METHODS);
if isempty(tag)
    fclose(fid);
    error(me,'Modalities not found');
end
inv.methods = tag.data;
%
tag = find_tag(invs,FIFF.FIFF_MNE_SOURCE_ORIENTATION);
if isempty(tag)
    fclose(fid);
    error(me,'Source orientation constraints not found');
end
inv.source_ori = tag.data;
%
tag = find_tag(invs,FIFF.FIFF_MNE_SOURCE_SPACE_NPOINTS);
if isempty(tag)
    fclose(fid);
    error(me,'Number of sources not found');
end
inv.nsource = tag.data;
inv.nchan   = 0;
%
%   Coordinate frame
%
tag = find_tag(invs,FIFF.FIFF_MNE_COORD_FRAME);
if isempty(tag)
    fclose(fid);
    error(me,'Coordinate frame tag not found');
end
inv.coord_frame = tag.data;
%
%   The actual source orientation vectors
%
tag = find_tag(invs,FIFF.FIFF_MNE_INVERSE_SOURCE_ORIENTATIONS);
if isempty(tag)
    fclose(fid);
    error(me,'Source orientation information not found');
end
inv.source_nn   = tag.data;
fprintf(1,'[done]\n');
%
%   The SVD decomposition...
%
fprintf(1,'\tReading inverse operator decomposition...');
tag = find_tag(invs,FIFF.FIFF_MNE_INVERSE_SING);
if isempty(tag)
    fclose(fid);
    error(me,'Singular values not found');
end
inv.sing  = tag.data;
inv.nchan = length(inv.sing);
%
%   The eigenleads and eigenfields
%
inv.eigen_leads_weighted = false;
try
   inv.eigen_leads = fiff_read_named_matrix(fid,invs,FIFF.FIFF_MNE_INVERSE_LEADS);
catch
   inv.eigen_leads_weighted = true;
   try
      inv.eigen_leads = fiff_read_named_matrix(fid,invs,FIFF.FIFF_MNE_INVERSE_LEADS_WEIGHTED);
   catch
      error(me,'%s',mne_omit_first_line(lasterr));
   end
end
%
%   Having the eigenleads as columns is better for the inverse calculations
%
inv.eigen_leads = mne_transpose_named_matrix(inv.eigen_leads);
try
    inv.eigen_fields = fiff_read_named_matrix(fid,invs,FIFF.FIFF_MNE_INVERSE_FIELDS);
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end
fprintf(1,'[done]\n');
%
%   Read the covariance matrices
%
try 
    inv.noise_cov = mne_read_cov(fid,invs,FIFF.FIFFV_MNE_NOISE_COV);
    fprintf('\tNoise covariance matrix read.\n');
catch
    fclose(fid);
    error(me,'%s',mne_omit_first_line(lasterr));
end
try 
    inv.source_cov = mne_read_cov(fid,invs,FIFF.FIFFV_MNE_SOURCE_COV);
    fprintf('\tSource covariance matrix read.\n');
catch
    fclose(fid);
    error(me,'%s',mne_omit_first_line(lasterr));
end
%
%   Read the various priors
%
try 
    inv.orient_prior = mne_read_cov(fid,invs,FIFF.FIFFV_MNE_ORIENT_PRIOR_COV);
    fprintf('\tOrientation priors read.\n');
catch
    inv.orient_prior = [];
end
try 
    inv.depth_prior = mne_read_cov(fid,invs,FIFF.FIFFV_MNE_DEPTH_PRIOR_COV);
    fprintf('\tDepth priors read.\n');
catch
    inv.depth_prior = [];
end
try 
    inv.fmri_prior = mne_read_cov(fid,invs,FIFF.FIFFV_MNE_FMRI_PRIOR_COV);
    fprintf('\tfMRI priors read.\n');
catch
    inv.fmri_prior = [];
end
%
%   Read the source spaces
%
try
    inv.src = mne_read_source_spaces(fid,false,tree);
catch
    fclose(fid);
    error(me,'Could not read the source spaces (%s)',mne_omit_first_line(lasterr));
end
for k = 1:length(inv.src)
   inv.src(k).id = mne_find_source_space_hemi(inv.src(k));
end
%
%   Get the MRI <-> head coordinate transformation
%
tag = find_tag(parent_mri,FIFF.FIFF_COORD_TRANS);
if isempty(tag)
    fclose(fid);
    error(me,'MRI/head coordinate transformation not found');
else
    mri_head_t = tag.data;
    if mri_head_t.from ~= FIFF.FIFFV_COORD_MRI || mri_head_t.to ~= FIFF.FIFFV_COORD_HEAD
        mri_head_t = fiff_invert_transform(mri_head_t);
        if mri_head_t.from ~= FIFF.FIFFV_COORD_MRI || mri_head_t.to ~= FIFF.FIFFV_COORD_HEAD
            fclose(fid);
            error(me,'MRI/head coordinate transformation not found');
        end
    end
end
inv.mri_head_t  = mri_head_t;
%
%   Transform the source spaces to the correct coordinate frame
%   if necessary
%
if inv.coord_frame ~= FIFF.FIFFV_COORD_MRI && ...
        inv.coord_frame ~= FIFF.FIFFV_COORD_HEAD
    fclose(fid);
    error(me,'Only inverse solutions computed in MRI or head coordinates are acceptable');
end
%
%  Number of averages is initially one
%
inv.nave = 1;
%
%  We also need the SSP operator
%
inv.projs     = fiff_read_proj(fid,tree);
%
%  Some empty fields to be filled in later
%
inv.proj      = [];      %   This is the projector to apply to the data
inv.whitener  = [];      %   This whitens the data
inv.reginv    = [];      %   This the diagonal matrix implementing
                         %   regularization and the inverse
inv.noisenorm = [];      %   These are the noise-normalization factors
%
nuse = 0;
for k = 1:length(inv.src)
   try
      inv.src(k) = mne_transform_source_space_to(inv.src(k),inv.coord_frame,mri_head_t);
   catch
      fclose(fid);
      error(me,'Could not transform source space (%s)',mne_omit_first_line(lasterr));
   end
   nuse = nuse + inv.src(k).nuse;
end
fprintf(1,'\tSource spaces transformed to the inverse solution coordinate frame\n');
%
%   Done!
%
fclose(fid);

return;

    function [tag] = find_tag(node,findkind)

        for p = 1:node.nent
            if node.dir(p).kind == findkind
                tag = fiff_read_tag(fid,node.dir(p).pos);
                return;
            end
        end
        tag = [];
        return;
    end
        
end

 
 
