function mne_write_cov(fid,cov)
%
%
%   mne_write_cov(fid,cov)
%
%   Write a covariance matrix to an open file
%
%   fid     - an open file id
%   cov     - the covariance matrix to write
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.4  2007/12/10 01:02:55  msh
%   Fixed writing of a diagonal covariance matrix
%
%   Revision 1.3  2006/05/05 03:50:40  msh
%   Added routines to compute L2-norm inverse solutions.
%   Added mne_write_inverse_sol_stc to write them in stc files
%   Several bug fixes in other files
%
%   Revision 1.2  2006/05/03 18:53:06  msh
%   Approaching Matlab 6.5 backward compatibility
%
%   Revision 1.1  2006/04/29 12:44:10  msh
%   Added covariance matrix writing routines.
%
%

me='MNE:mne_write_cov';

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end

fiff_start_block(fid,FIFF.FIFFB_MNE_COV);
%
%   Dimensions etc.
%
fiff_write_int(fid,FIFF.FIFF_MNE_COV_KIND,cov.kind);
fiff_write_int(fid,FIFF.FIFF_MNE_COV_DIM,cov.dim);
if cov.nfree > 0
    fiff_write_int(fid,FIFF.FIFF_MNE_COV_NFREE,cov.nfree);
end
%
%   Channel names
%
if ~isempty(cov.names)
    fiff_write_name_list(fid,FIFF.FIFF_MNE_ROW_NAMES,cov.names);
end
%
%   Data
%
if cov.diag
   fiff_write_double(fid,FIFF.FIFF_MNE_COV_DIAG,cov.data);
else
   if issparse(cov.data)
      fiff_write_float_sparse_rcs(fid,FIFF.FIFF_MNE_COV,cov.data);
   else
      q = 1;
      vals = zeros(cov.dim*(cov.dim+1)/2,1);
      for j = 1:cov.dim
          for k = 1:j
              vals(q) = cov.data(j,k);
              q = q + 1;
          end
      end
      fiff_write_double(fid,FIFF.FIFF_MNE_COV,vals);
   end
end
%
%   Eigenvalues and vectors if present
%
if ~isempty(cov.eig) && ~isempty(cov.eigvec)
    fiff_write_float_matrix(fid,FIFF.FIFF_MNE_COV_EIGENVECTORS,cov.eigvec);
    fiff_write_double(fid,FIFF.FIFF_MNE_COV_EIGENVALUES,cov.eig);
end
%
%   Projection operator
%
fiff_write_proj(fid,cov.projs);
%
%   Bad channels
%
if ~isempty(cov.bads)
    fiff_start_block(fid,FIFF.FIFFB_MNE_BAD_CHANNELS);
    fiff_write_name_list(fid,FIFF.FIFF_MNE_CH_NAME_LIST,cov.bads);
    fiff_end_block(fid,FIFF.FIFFB_MNE_BAD_CHANNELS);
end
%
%   Done!
%
fiff_end_block(fid,FIFF.FIFFB_MNE_COV);
