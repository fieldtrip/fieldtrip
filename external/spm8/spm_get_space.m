function M = spm_get_space(P,M)
% Get/set the voxel-to-world mapping of an image
% FORMAT M = spm_get_space(P)
%            spm_get_space(P,M)
% M - voxel-to-world mapping
% P - image filename
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_get_space.m 4208 2011-02-22 12:04:18Z guillaume $


[pth,nam,ext,num] = spm_fileparts(P);
if ~isempty(num), n = str2num(num(2:end)); else n = [1 1]; end
P = fullfile(pth,[nam ext]);

N = nifti(P);
if nargin==2
    N.mat_intent = 'Aligned';
    if n(1)==1

        % Ensure volumes 2..N have the original matrix
        if size(N.dat,4)>1 && sum(sum((N.mat-M).^2))>1e-8
            M0 = N.mat;
            if ~isfield(N.extras,'mat')
                N.extras.mat = zeros([4 4 size(N.dat,4)]);
            else
                if size(N.extras.mat,4)<size(N.dat,4)
                    N.extras.mat(:,:,size(N.dat,4)) = zeros(4);
                end
            end
            for i=2:size(N.dat,4)
                if sum(sum(N.extras.mat(:,:,i).^2))==0,
                    N.extras.mat(:,:,i) = M0;
                end
            end
        end

        N.mat = M;
        if strcmp(N.mat0_intent,'Aligned'), N.mat0 = M; end
        if ~isempty(N.extras) && isstruct(N.extras) && isfield(N.extras,'mat') &&...
            size(N.extras.mat,3)>=1
            N.extras.mat(:,:,n(1)) = M;
        end
    else
        N.extras.mat(:,:,n(1)) = M;
    end
    create(N);
else
    if ~isempty(N.extras) && isstruct(N.extras) && isfield(N.extras,'mat') &&...
        size(N.extras.mat,3)>=n(1) && sum(sum(N.extras.mat(:,:,n(1)).^2))
        M = N.extras.mat(:,:,n(1));
    else
        M = N.mat;
    end
end
