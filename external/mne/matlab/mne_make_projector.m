function [proj,nproj,U] = mne_make_projector(projs,ch_names,bads)
%
% [proj,nproj,U] = mne_make_projector(projs,ch_names,bads)
%
% proj     - The projection operator to apply to the data
% nproj    - How many items in the projector
% U        - The orthogonal basis of the projection vectors (optional)
%
% Make an SSP operator
%
% projs    - A set of projection vectors
% ch_names - A cell array of channel names
% bads     - Bad channels to exclude
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.11  2008/11/22 17:06:42  msh
%   Optionally return U
%
%   Revision 1.10  2006/05/05 19:37:47  msh
%   Fixed error in mne_make_projector.
%   Better detection of small eigenvalues for the projector.
%
%   Revision 1.9  2006/05/05 03:50:40  msh
%   Added routines to compute L2-norm inverse solutions.
%   Added mne_write_inverse_sol_stc to write them in stc files
%   Several bug fixes in other files
%
%   Revision 1.8  2006/04/26 00:43:22  msh
%   Fixed errors in mne_make_projector related to vectors which do not affect the data
%
%   Revision 1.7  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.6  2006/04/21 14:23:16  msh
%   Further improvements in raw data reading
%
%   Revision 1.5  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.4  2006/04/17 15:01:34  msh
%   More small improvements.
%
%   Revision 1.3  2006/04/15 12:21:00  msh
%   Several small improvements
%
%   Revision 1.2  2006/04/14 15:49:49  msh
%   Improved the channel selection code and added ch_names to measurement info.
%
%   Revision 1.1  2006/04/14 03:30:49  msh
%   Added mne_make_projector
%
%

me='MNE:mne_make_projector';

if nargin == 2
    bads = [];
elseif nargin ~= 3
    error(me,'Incorrect number of arguments');
end

nchan = length(ch_names);
if nchan == 0
    error(me,'No channel names specified');
end

proj  = eye(nchan,nchan);
nproj = 0;
U     = [];
%
%   Check trivial cases first
%
if isempty(projs)
    return;
end

nactive = 0;
nvec    = 0;
for k = 1:length(projs)
    if projs(k).active
        nactive = nactive + 1;
        nvec = nvec + projs(k).data.nrow;
    end
end

if nactive == 0
    return;
end
%
%   Pick the appropriate entries
%
vecs = zeros(nchan,nvec);
nvec = 0;
nonzero = 0;
for k = 1:length(projs)
    if projs(k).active
        one = projs(k);
        sel = [];
        vecsel = [];
        if length(one.data.col_names) ~= length(unique(one.data.col_names))
            error(me,'Channel name list in projection item %d contains duplicate items',k);
        end
        %
        % Get the two selection vectors to pick correct elements from
        % the projection vectors omitting bad channels
        %
        p = 0;
        for c = 1:nchan
            match = strmatch(ch_names{c},one.data.col_names);
            if ~isempty(match) && isempty(strmatch(ch_names{c},bads))
                p = p + 1;
                sel(p)    = c;
                vecsel(p) = match(1);
            end
        end
        %
        % If there is something to pick, pickit
        %
        if ~isempty(sel)
            for v = 1:one.data.nrow
                vecs(sel,nvec+v) = one.data.data(v,vecsel)';
            end
            %
            %   Rescale for more straightforward detection of small singular values
            %
            for v = 1:one.data.nrow
                onesize = sqrt(vecs(:,nvec+v)'*vecs(:,nvec+v));
                if onesize > 0
                    vecs(:,nvec+v) = vecs(:,nvec+v)/onesize;
                    nonzero = nonzero + 1;
                end
            end
            nvec = nvec + one.data.nrow;
        end
    end
end
%
%   Check whether all of the vectors are exactly zero
%
if nonzero == 0
    return;
end
%
%   Reorthogonalize the vectors
%
[U,S,V] = svd(vecs(:,1:nvec),0);
S = diag(S);
%
%   Throw away the linearly dependent guys
%
for k = 1:nvec
    if S(k)/S(1) < 1e-2
        nvec = k;
        break;
    end
end
U = U(:,1:nvec);
%
%   Here is the celebrated result
%
proj  = proj - U*U';
nproj = nvec;

return;

end
