function [cov] = mne_read_cov(fid,node,cov_kind)
%
% [cov] = mne_read_cov(fid,node,kind)
%
% Reads a covariance matrix from a fiff file
%
% fid       - an open file descriptor
% node      - look for the matrix in here
% cov_kind  - what kind of a covariance matrix do we want?
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.6  2008/09/26 20:49:26  msh
%   Enabled reading of a sparse covariance matrix
%
%   Revision 1.5  2007/01/29 21:21:22  msh
%   Added reading of the additional source prior covariances.
%
%   Revision 1.4  2006/05/03 18:53:06  msh
%   Approaching Matlab 6.5 backward compatibility
%
%   Revision 1.3  2006/04/29 12:44:10  msh
%   Added covariance matrix writing routines.
%
%   Revision 1.2  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.1  2006/04/20 21:50:04  msh
%   Added mne_read_cov.m
%
%

me='MNE:mne_read_cov';

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end
%
%   Find all covariance matrices
%
covs = fiff_dir_tree_find(node,FIFF.FIFFB_MNE_COV);
if isempty(covs)
    error(me,'No covariance matrices found');
end
%
%   Is any of the covariance matrices a noise covariance
%
for p = 1:length(covs)
    tag = find_tag(covs(p),FIFF.FIFF_MNE_COV_KIND);
    if ~isempty(tag) && tag.data == cov_kind
        this = covs(p);
        %
        %   Find all the necessary data
        %
        tag = find_tag(this,FIFF.FIFF_MNE_COV_DIM);
        if isempty(tag)
            error(me,'Covariance matrix dimension not found');
        end
        dim = tag.data;
        tag = find_tag(this,FIFF.FIFF_MNE_COV_NFREE);
        if isempty(tag)
            nfree = -1;
        else
            nfree = tag.data;
        end
        tag = find_tag(this,FIFF.FIFF_MNE_ROW_NAMES);
        if isempty(tag)
            names = [];
        else
            names = fiff_split_name_list(tag.data);
            if size(names,2) ~= dim
                error(me,'Number of names does not match covariance matrix dimension');
            end
        end
        tag = find_tag(this,FIFF.FIFF_MNE_COV);
        if isempty(tag)
            tag = find_tag(this,FIFF.FIFF_MNE_COV_DIAG);
            if isempty(tag)
                error(me,'No covariance matrix data found');
            else
                %
                %   Diagonal is stored
                %
                data = tag.data;
                diagmat = true;
                fprintf('\t%d x %d diagonal covariance (kind = %d) found.\n',dim,dim,cov_kind);
            end
        else
            if ~issparse(tag.data)
                %
                %   Lower diagonal is stored
                %
                vals = tag.data;
                data = zeros(dim,dim);
                % XXX : should remove for loops
                q = 1;
                for j = 1:dim
                    for k = 1:j
                        data(j,k) = vals(q);
                        q = q + 1;
                    end
                end
                for j = 1:dim
                    for k = j+1:dim
                        data(j,k) = data(k,j);
                    end
                end
                diagmat = false;
                fprintf('\t%d x %d full covariance (kind = %d) found.\n',dim,dim,cov_kind);
            else
                diagmat = false;
                data = tag.data;
                fprintf('\t%d x %d sparse covariance (kind = %d) found.\n',dim,dim,cov_kind);
            end
        end
        %
        %   Read the possibly precomputed decomposition
        %
        tag1 = find_tag(this,FIFF.FIFF_MNE_COV_EIGENVALUES);
        tag2 = find_tag(this,FIFF.FIFF_MNE_COV_EIGENVECTORS);
        if ~isempty(tag1) && ~isempty(tag2)
            eig    = tag1.data;
            eigvec = tag2.data;
        else
            eig    = [];
            eigvec = [];
        end
        %
        %   Read the projection operator
        %
        projs = fiff_read_proj(fid,this);
        %
        %   Read the bad channel list
        %
        bads = fiff_read_bad_channels(fid,this);
        %
        %   Put it together
        %
        cov.kind   = cov_kind;
        cov.diag   = diagmat;
        cov.dim    = dim;
        cov.names  = names;
        cov.data   = data;
        cov.projs  = projs;
        cov.bads   = bads;
        cov.nfree  = nfree;
        cov.eig    = eig;
        cov.eigvec = eigvec;
        %
        return;
    end
end

error(me,'Did not find the desired covariance matrix');

return;

    function [tag] = find_tag(node,findkind)
        
        for pp = 1:node.nent
            kind = node.dir(pp).kind;
            pos  = node.dir(pp).pos;
            if kind == findkind
                tag = fiff_read_tag(fid,pos);
                return;
            end
        end
        tag = [];
        return
    end

end
