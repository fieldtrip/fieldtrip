function [res] = mne_ex_compute_inverse(fname_data,setno,fname_inv,nave,lambda2,dSPM,sLORETA)
%
% [res] = mne_ex_compute_inverse(fname_data,setno,fname_inv,nave,lambda2,dSPM,sLORETA)
%
% An example on how to compute a L2-norm inverse solution
% Actual code using these principles might be different because 
% the inverse operator is often reused across data sets.
%
%
% fname_data  - Name of the data file
% setno       - Data set number
% fname_inv   - Inverse operator file name
% nave        - Number of averages (scales the noise covariance)
%               If negative, the number of averages in the data will be
%               used
% lambda2     - The regularization factor
% dSPM        - do dSPM?
% sLORETA     - do sLORETA?
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.2  2008/05/26 10:49:26  msh
%   Update to incorporate the already weighted lead field basis
%
%   Revision 1.1  2006/05/05 03:50:40  msh
%   Added routines to compute L2-norm inverse solutions.
%   Added mne_write_inverse_sol_stc to write them in stc files
%   Several bug fixes in other files
%
%

me='MNE:mne_ex_compute_inverse';

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end

if nargin ~= 6 && nargin ~= 7
   error(me,'Incorrect number of arguments'); 
end

if nargin == 6
   sLORETA = false;
end
%
%   Read the data first
%
data = fiff_read_evoked(fname_data,setno);
%
%   Then the inverse operator
%
inv = mne_read_inverse_operator(fname_inv);
%
%   Set up the inverse according to the parameters
%
if nave < 0
    nave = data.evoked.nave;
end
inv = mne_prepare_inverse_operator(inv,nave,lambda2,dSPM,sLORETA);
%
%   Pick the correct channels from the data
%
data = fiff_pick_channels_evoked(data,inv.noise_cov.names);
fprintf(1,'Picked %d channels from the data\n',data.info.nchan);
fprintf(1,'Computing inverse...');
%
%   Simple matrix multiplication followed by combination of the 
%   three current components
%
%   This does all the data transformations to compute the weights for the
%   eigenleads
%   
trans = diag(sparse(inv.reginv))*inv.eigen_fields.data*inv.whitener*inv.proj*double(data.evoked(1).epochs);
%
%   Transformation into current distributions by weighting the eigenleads
%   with the weights computed above
%
if inv.eigen_leads_weighted
   %
   %     R^0.5 has been already factored in
   %
   fprintf(1,'(eigenleads already weighted)...');
   sol   = inv.eigen_leads.data*trans;
else
   %
   %     R^0.5 has to factored in
   %
   fprintf(1,'(eigenleads need to be weighted)...');
   sol   = diag(sparse(sqrt(inv.source_cov.data)))*inv.eigen_leads.data*trans;
end
   
if inv.source_ori == FIFF.FIFFV_MNE_FREE_ORI
    fprintf(1,'combining the current components...');
    sol1 = zeros(size(sol,1)/3,size(sol,2));
    for k = 1:size(sol,2)
        sol1(:,k) = sqrt(mne_combine_xyz(sol(:,k)));
    end
    sol = sol1;
end
if dSPM
    fprintf(1,'(dSPM)...');
    sol = inv.noisenorm*sol;
elseif sLORETA
    fprintf(1,'(sLORETA)...');
    sol = inv.noisenorm*sol;
end
res.inv   = inv;
res.sol   = sol;
res.tmin  = double(data.evoked(1).first)/data.info.sfreq;
res.tstep = 1/data.info.sfreq;
fprintf(1,'[done]\n');

return;
end

