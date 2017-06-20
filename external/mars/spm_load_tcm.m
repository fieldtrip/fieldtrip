function tcm = spm_load_tcm(C)
% tcm = spm_load_tcm(C)
%
% Adapted from SPM function spm_load_priors(V)
% Load the local version of tissue correlation maps for segmentation,
% mainly getting the B-spline coefficients for subsequent sampling using
% B-spline interpolation.
%
% C: local TCM matrix (dim1 x dim2 x dim3 x K x K x 6, dim1~dim3 are same
% as the image size, K is number of tissue types, 6 means 6-neighborhood)
% tcm: a structure for TCM
%
% This function is intended to be used in conjunction with spm_sample_tcm.
% tcm = spm_load_tcm(C);
% C_slice = spm_sample_tcm(tcm,X,Y,Z);
%____________________________________________________________________________
%
% John Ashburner
% Yu (Andy) Huang, 2013-06-04

tiny = eps*eps; % 1e-3; % ANDY 2013-06-05

deg = 1;

% if ~isstruct(V), V  = spm_vol(V); end;
% spm_check_orientations(V);

% =========================NOTE!!!!!!=====================================
% TCM is saved as .mat file only (NOT as NIFTI file).
% Because if NIFTI, it'll have huge size. The header info is NOT needed as
% long as TCM is generated from the same dataset as TPM (then it's already in correct orientation)
% If the dataset that was used to generate TCM does NOT have same header
% info as the TPM NIFTI file, then this interpolation will NOT work
% properly for TCM.
% Therefore, TCM should be generated from the same dataset that was used to
% get TPM, at least header info should be same. % ANDY 2013-06-05
% =========================NOTE!!!!!!=====================================

% tpm.V = V;
% tpm.M = tpm.V(1).mat;

d = size(C);
% Kb = numel(tpm.V);
tcm.dat = cell(d(4:6));

% spm_progress_bar('Init',tpm.V(1).dim(3),'Loading priors','Planes loaded');
% for i=1:tpm.V(1).dim(3)
%     M         = spm_matrix([0 0 i]);
%     s         = zeros(tpm.V(1).dim(1:2));
%     for k1=1:Kb
%         tmp                = spm_slice_vol(tpm.V(k1),M,tpm.V(1).dim(1:2),0);
%         tpm.dat{k1}(:,:,i) = max(min(tmp,1),0);
%         s                  = s + tmp;
%     end;
%     t = s>1;
%     if any(t)
%         for k1=1:Kb
%             tmp           = tpm.dat{k1}(:,:,i);
%             tmp(t)        = tmp(t)./s(t); % normalize the loaded TPM
%             % eTPM.nii will be normalized here. Others (bTPM, cTPM, cTPMthresh) were normalized before loading.
%             tpm.dat{k1}(:,:,i) = tmp;
%         end;
%     end;
%     spm_progress_bar('Set',i);
% end;

tcm.bg1 = zeros(d(4:6));
tcm.bg2 = zeros(d(4:6));
for k=1:d(6)
    for j=1:d(5)
        for i=1:d(4)
            tcm.bg1(i,j,k) = mean(mean(C(:,:,1,i,j,k)));
            tcm.bg2(i,j,k) = mean(mean(C(:,:,end,i,j,k)));
            tcm.dat{i,j,k} = spm_bsplinc(log(C(:,:,:,i,j,k)+tiny),[deg deg deg  0 0 0]);
            % get B-spline coefficients for subsequent sampling using B-spline interpolation

%             tcm.dat{i,j,k} = spm_bsplinc(log(C(2:(end-1),2:(end-1),2:(end-1),i,j,k)+tiny),[deg deg deg  0 0 0]);
%             % here use the 2nd slice and the penultimate slice, because for TCM, it has edge effects, on the edges, some neighbors are not defined.
%             % therefore, edge slices are all discarded when computing B-spline coefficients % ANDY 2013-06-05
        end
    end
end

% tcm.tiny = tiny;
tcm.deg  = deg+1;
% spm_progress_bar('Clear');
return;