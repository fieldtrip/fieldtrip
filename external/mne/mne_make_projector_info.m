function [proj,nproj] = mne_make_projector_info(info)
%
% [proj,nproj] = mne_make_projector_info(info)
%
% Make an SSP operator using the meas info
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.1  2006/05/05 10:52:33  msh
%   Added missing file.
%
%

me='MNE:mne_make_projector_info';

if nargin ~= 1
   error(me,'Incorrect number of arguments');
end

[ proj, nproj ] = mne_make_projector(info.projs,info.ch_names,info.bads);

return;

end
