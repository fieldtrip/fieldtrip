function spm_dartel_dotprods(job)
% Generate a kernel from dot-products of images
% FORMAT spm_dartel_dotprods(job)
% job.images  - Images to use
% job.dotprod - Part of filename for results
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dartel_dotprods.m 4492 2011-09-16 12:11:09Z guillaume $

P      = strvcat(job.images);
[pth,nam,ext] = fileparts(job.dotprod);
ofname = fullfile(pwd,['dp_' nam '.mat']);

N = nifti(P);
n = numel(N);
dm= size(N(1).dat);
dat=cell(1,numel(N));
for i=1:numel(N),
    dat{i} = reshape(N(i).dat,[prod(dm),1]);
end
Phi = zeros(n,n);

if isfield(job,'weight') && ~isempty(job.weight),
    Pmsk = strvcat(job.weight);
    Nmsk = nifti(Pmsk);
    msk  = Nmsk.dat;
    dmsk = size(msk);
    if any(dmsk(1:3) ~= dm(1:3)),
        error('Wrong sized weighting image.');
    end
    msk = reshape(msk,[prod(dmsk),1]);
    if numel(dmsk)==3,
        msk1 = msk;
        for i=2:prod(dm(4:end)),
            msk = [msk;msk1];
        end
    end
end

mem = 32*1024*1024;  % Mbytes of RAM to use
bs  = ceil(mem/8/n); % Block size
nd  = prod(dm);
nblock = ceil(prod(dm)/bs);
spm_progress_bar('Init',nblock,...
                 'Generating kernel','Blocks complete');
for k=1:nblock,
    o = bs*(k-1)+(1:bs);
    o = o(o<nd);
    if exist('msk','var'),
        wt  = msk(o);
        tmp = wt>0;
        o   = o(tmp);
        wt  = wt(tmp);
    end
    if ~isempty(o),
        X = zeros(numel(o),numel(dat));
        for i=1:n,
            tmp    = dat{i}(o);
            tmp(~isfinite(tmp)) = 0;
            if exist('wt','var'), tmp = tmp.*wt; end
            X(:,i) = tmp;
        end
        Phi    = Phi + X'*X;
        clear X
    end
    spm_progress_bar('Set',k);
end
spm_progress_bar('Clear');
input = job;
typ   = 'images';
save(ofname,'Phi','input','typ', spm_get_defaults('mat.format'));

