function [inv] = mne_prepare_inverse_operator(orig,nave,lambda2,dSPM)
%
% [inv] = mne_prepare_inverse_operator(orig,nave,lambda2,dSPM)
%
% Prepare for actually computing the inverse
%
% orig        - The inverse operator structure read from a file
% nave        - Number of averages (scales the noise covariance)
% lambda2     - The regularization factor
% dSPM        - Compute the noise-normalization factors for dSPM?
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.4  2009/03/09 11:22:23  msh
%   Fixed the noise-normalization factors for the case where eigen leads come
%   already weighted
%
%   Revision 1.3  2008/05/26 10:49:26  msh
%   Update to incorporate the already weighted lead field basis
%
%   Revision 1.2  2006/05/05 19:37:47  msh
%   Fixed error in mne_make_projector.
%   Better detection of small eigenvalues for the projector.
%
%   Revision 1.1  2006/05/05 03:50:40  msh
%   Added routines to compute L2-norm inverse solutions.
%   Added mne_write_inverse_sol_stc to write them in stc files
%   Several bug fixes in other files
%
%

me='MNE:mne_prepare_inverse_operator';

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

if nargin ~= 4
    error(me,'Wrong number of arguments');
end

if nave <= 0
    error(me,'The number of averages should be positive');
end
fprintf(1,'Preparing the inverse operator for use...\n');
inv = orig;
%
%   Scale some of the stuff
%
scale               = double(inv.nave)/double(nave);
inv.noise_cov.data  = scale*inv.noise_cov.data;
inv.noise_cov.eig   = scale*inv.noise_cov.eig;
inv.source_cov.data = scale*inv.source_cov.data;
%
if inv.eigen_leads_weighted
    inv.eigen_leads.data = sqrt(scale)*inv.eigen_leads.data;
end
%
fprintf(1,'\tScaled noise and source covariance from nave = %d to nave = %d\n',inv.nave,nave);
inv.nave = nave;
%
%   Create the diagonal matrix for computing the regularized inverse
%
inv.reginv = inv.sing./(inv.sing.*inv.sing + lambda2);
fprintf(1,'\tCreated the regularized inverter\n');
%
%   Create the projection operator
%
[ inv.proj, ncomp ] = mne_make_projector(inv.projs,inv.noise_cov.names);
if ncomp > 0
    
    fprintf(1,'\tCreated an SSP operator (subspace dimension = %d)\n',ncomp);
end
%
%   Create the whitener
%
inv.whitener = zeros(inv.noise_cov.dim);
if inv.noise_cov.diag == 0
    %
    %   Omit the zeroes due to projection
    %
    nnzero = 0;
    for k = ncomp+1:inv.noise_cov.dim
        if inv.noise_cov.eig(k) > 0
            inv.whitener(k,k) = 1.0/sqrt(inv.noise_cov.eig(k));
            nnzero = nnzero + 1;
        end
    end
    %
    %   Rows of eigvec are the eigenvectors
    %
    inv.whitener = inv.whitener*inv.noise_cov.eigvec;
    fprintf(1,'\tCreated the whitener using a full noise covariance matrix (%d small eigenvalues omitted)\n', ...
        inv.noise_cov.dim - nnzero);
else
    %
    %   No need to omit the zeroes due to projection
    %
    for k = 1:inv.noise_cov.dim
        inv.whitener(k,k) = 1.0/sqrt(inv.noise_cov.data(k));
    end
    fprintf(1,'\tCreated the whitener using a diagonal noise covariance matrix (%d small eigenvalues discarded)\n',ncomp);
end
%
%   Finally, compute the noise-normalization factors
%
if dSPM
    fprintf(1,'\tComputing noise-normalization factors...');
    noise_norm = zeros(inv.eigen_leads.nrow,1);
    if inv.eigen_leads_weighted
        for k = 1:inv.eigen_leads.nrow
            one = inv.eigen_leads.data(k,:).*inv.reginv';
            noise_norm(k) = sqrt(one*one');
        end
    else
        for k = 1:inv.eigen_leads.nrow
            one = sqrt(inv.source_cov.data(k))*(inv.eigen_leads.data(k,:).*inv.reginv');
            noise_norm(k) = sqrt(one*one');
        end
    end
    %
    %   Compute the final result
    %
    if inv.source_ori == FIFF.FIFFV_MNE_FREE_ORI
        %
        %   The three-component case is a little bit more involved
        %   The variances at three consequtive entries must be squeared and
        %   added together
        %
        %   Even in this case return only one noise-normalization factor
        %   per source location
        %
        noise_norm = sqrt(mne_combine_xyz(noise_norm));
        %
        %   This would replicate the same value on three consequtive
        %   entries
        %
        %   noise_norm = kron(sqrt(mne_combine_xyz(noise_norm)),ones(3,1));
    end
    inv.noisenorm = diag(sparse(1./abs(noise_norm)));
    fprintf(1,'[done]\n');
else
    inv.noisenorm = [];
end

return;

end
