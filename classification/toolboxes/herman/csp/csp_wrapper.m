
function [W,d] = csp_wrapper(filtered_data1,filtered_data2,k_pat,mode)
%
% CSP_WRAPPER is an auxiliary function used by csp_train (on its own or as
% part of CSPPROCESSOR object). It extracts CSP filter/patterns from two-class data.
% Essentially, it is a wrapper for biosig version of CSP (see csp.m).
%
% Use as
%         1) Used by CSP_TRAIN (csp_train.m)
%         2) As part of CSPPROCESSOR object (see cspprocessor.m) - RECOMMENDED
%         3) [W,d] = csp_wrapper(filtered_data1,filtered_data2,k_pat,mode)
%
% INPUT
%           filtered_data1, 
%           filtered_data2 - three-dimensional arrays [trial x chan x time] 
%                            of data corresponding to two classes
%           k_pat          - number of CSP filters (2*k_pat)
%           mode           - 'CSP0' or 'CSP3' (see csp.m or csp_train.m)
%
% REMARKS
%          1) Data can be prefiltered in the desirable frequency band (see PREPROCESSING)
%          2) Covariance is calculated along the 3rd dimension and averaged over the 1st one
%          3) Requires ecovm.m and scp.m from BIOSIG
%
% OUTPUT 
%           W - arrays of CSP filters (2*k_pat column vectors) derived from class1 and class2
%           d - corresponding eigenvalues
%
% SEE ALSO
%           csp.m
%           csp_test.m
%           csp_train.m
%           cspprocessor.m
%

% Pawel Herman, 2009


if nargin<3, k_pat= 1; end
if nargin<4, mode = 'CSP0'; end

ECOV1= 0; ECOV2 = 0;
for i=1:size(filtered_data1,1), ECOV1 = ECOV1 + ecovm(reshape(filtered_data1(i,:,:),size(filtered_data1,2),size(filtered_data1,3))'); end
for i=1:size(filtered_data2,1), ECOV2 = ECOV2 + ecovm(reshape(filtered_data2(i,:,:),size(filtered_data2,2),size(filtered_data2,3))'); end
eCOV(1,:,:) = ECOV1;   
eCOV(2,:,:) = ECOV2;

[W d] = csp(eCOV,k_pat,mode);



