function b0 = spm_load_priors(B)
% Loads the tissue probability maps for segmentation
% FORMAT b0 = spm_load_priors(B)
% B  - structures of image volume information (or filenames)
% b0 - a cell array of tissue probabilities
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_load_priors.m 4873 2012-08-30 19:06:26Z john $


% deg = 3;
lm = 0;
if ~isstruct(B), B  = spm_vol(B); end
Kb = length(B);
b0 = cell(Kb,1);
for k1=1:(Kb)
    b0{k1} = zeros(B(1).dim(1:3));
end

spm_progress_bar('Init',B(1).dim(3),'Loading priors','Planes loaded');
for i=1:B(1).dim(3)
    M         = spm_matrix([0 0 i]);
    s         = zeros(B(1).dim(1:2));
    for k1=1:Kb
        tmp           = spm_slice_vol(B(k1),M,B(1).dim(1:2),0)*(1-lm*2)+lm;
        b0{k1}(:,:,i) = max(min(tmp,1),0);
        s             = s + tmp;
    end
    t = s>1;
    if any(any(t))
        for k1=1:Kb
            tmp           = b0{k1}(:,:,i);
            tmp(t)        = tmp(t)./s(t);
            b0{k1}(:,:,i) = tmp;
        end
    end
    s(t) = 1;
    b0{Kb+1}(:,:,i) = max(min(1-s,1),0);
    spm_progress_bar('Set',i);
end
%for k1=1:Kb+1
%    b0{k1} = spm_bsplinc(log(b0{k1}),[deg deg deg  0 0 0]);
%end
spm_progress_bar('Clear');
