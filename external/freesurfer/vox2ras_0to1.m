function M1 = vox2ras_0to1(M0)
% M1 = vox2ras_0to1(M0)
%
% Converts a 0-based vox2ras matrix to 1-based, ie,
% Pxyz = M0*[c r s 1]' = M1*[(c+1) (r+1) (s+1) 1]'
%


%
% vox2ras_0to1.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:13 $
%    $Revision: 1.3 $
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


M1 = [];

if(nargin ~= 1)
  fprintf('M1 = vox2ras_0to1(M0)\n');
  return;
end

Q = zeros(4);
Q(1:3,4) = ones(3,1);

M1 = inv(inv(M0)+Q);

return;







