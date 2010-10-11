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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:10 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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







