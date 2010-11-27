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
%    $Date: 2007/01/10 22:55:10 $
%    $Revision$
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


M1 = [];

if(nargin ~= 1)
  fprintf('M1 = vox2ras_0to1(M0)\n');
  return;
end

Q = zeros(4);
Q(1:3,4) = ones(3,1);

M1 = inv(inv(M0)+Q);

return;







