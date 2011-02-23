function [cov] = mne_pick_channels_cov(orig,include,exclude)
%
% [cov] = mne_pick_channels_cov(orig,include,exclude)
%
% Pick desired channels from a covariance matrix
%
% orig      - The original covariance matrix
% include   - Channels to include (if empty, include all available)
% exclude   - Channels to exclude (if empty, do not exclude any)
%
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.9  2006/11/21 12:53:49  msh
%   Fixed error in picking
%
%   Revision 1.8  2006/06/29 22:12:28  msh
%   Fixed errors in channel picking
%
%   Revision 1.7  2006/05/03 18:53:05  msh
%   Approaching Matlab 6.5 backward compatibility
%
%   Revision 1.6  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.5  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.4  2006/04/15 12:21:00  msh
%   Several small improvements
%
%   Revision 1.3  2006/04/14 15:49:49  msh
%   Improved the channel selection code and added ch_names to measurement info.
%
%   Revision 1.2  2006/04/14 00:45:42  msh
%   Added channel picking and fiff_invert_transform
%
%   Revision 1.1  2006/04/13 17:05:45  msh
%   Added reading of bad channels to fiff_read_meas_info.m
%   Added mne_pick_channels_cov.m
%
%

me='MNE:mne_pick_channels_cov';

if nargin == 1
    cov = orig;
    if isempty(cov.eig) || isempty(cov.eigvec)
        decompose_eigen; 
    end
    return;
elseif nargin == 2
    exclude = [];
elseif nargin ~= 3
    error(me,'Incorrect number of arguments');
end

if isempty(include) && isempty(exclude)
    cov = orig;
    if isempty(cov.eig) || isempty(cov.eigvec)
        decompose_eigen; 
    end
    return;
end

if isempty(orig.names) 
    error(me,'Cannot pick from a covariance matrix without channel names');
end

cov  = orig;
%
%   First do the channels to be included
%
sel = fiff_pick_channels(cov.names,include,exclude);
if isempty(sel)
   error(me,'Nothing remains after picking');
end
%
%   Select the desired stuff
%
if cov.diag
   cov.data = cov.data(sel);
else
   cov.data = cov.data(:,sel);
   cov.data = cov.data(sel,:);
end
for p = 1:size(sel,2)
   names{p} = cov.names{sel(p)};
end
cov.names = names;
cov.dim   = length(cov.names);
%
%   Eigenvalues and vectors are no longer valid
%
decompose_eigen;
%
return;

    function decompose_eigen
        if cov.diag 
            cov.eig    = cov.data;
            cov.eigvec = eye(cov.dim);
        else
            [ cov.eigvec, cov.eig ] = eig(cov.data);
            %
            %   We use the convention that rows of eigvec are the
            %   eigenvectors
            %
            cov.eigvec = cov.eigvec';
            cov.eig    = diag(cov.eig);
        end
    end

end
