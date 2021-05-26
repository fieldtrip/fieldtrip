function M0 = vox2ras_1to0(M1)
% M0 = vox2ras_1to0(M1)
%
% Converts a 1-based vox2ras matrix to 0-based, ie,
% Pxyz = M0*[c r s 1]' = M1*[(c+1) (r+1) (s+1) 1]'
%
% See also: vox2ras_0to1
%


%
% vox2ras_1to0.m
%
% Original Author: Doug Greve
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


M0 = [];

if(nargin ~= 1)
  fprintf('M0 = vox2ras_1to0(M1)');
  return;
end

Q = zeros(4);
Q(1:3,4) = -ones(3,1);

M0 = inv(inv(M1)+Q);

return;







