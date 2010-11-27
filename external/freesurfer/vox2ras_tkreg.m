function T = vox2ras_tkregmtx(voldim, voxres)
% T = vox2ras_tkregmtx(voldim, voxres)
%   voldim = [ncols  nrows  nslices ...]
%   volres = [colres rowres sliceres ...]
%
% Computes the 0-based vox2ras transform of a volume
% compatible with the registration matrix produced
% by tkregister. May work with MINC xfm.
%


%
% vox2ras_tkreg.m
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

T = [];

if(nargin ~= 2)
  fprintf('T = vox2ras_tkregmtx(voldim, voxres)\n');
  return;
end

T = zeros(4,4);
T(4,4) = 1;

T(1,1) = -voxres(1);
T(1,4) = voxres(1)*voldim(1)/2;

T(2,3) = voxres(3);
T(2,4) = -voxres(3)*voldim(3)/2;


T(3,2) = -voxres(2);
T(3,4) = voxres(2)*voldim(2)/2;


return;
