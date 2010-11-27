function [b,c,d,x,y,z,qfac] = vox2rasToQform(vox2ras)
% [a,b,c,x,y,z,qfac] = vox2rasToQform(vox2ras)
%   
% Converts a vox2ras matrix to NIFTI qform parameters.
% Note: the vox2ras should be 6 DOF. This code mostly just 
% follows CH's mriToNiftiQform() in mriio.c.
% 
%  hdr.pixdim(1) = qfac;
%  hdr.quatern_b = b;
%  hdr.quatern_c = c;
%  hdr.quatern_d = d;
%  hdr.qoffset_x = x;
%  hdr.qoffset_y = y;
%  hdr.qoffset_z = z;
%  hdr.qform_code = NIFTI_XFORM_SCANNER_ANAT=1;
%


%
% vox2rasToQform.m
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


a = [];
if(nargin ~= 1)
  fprintf('[a,b,c,x,y,z,qfac] = vox2rasToQform(vox2ras)\n');
  return;
end

x = vox2ras(1,4);
y = vox2ras(2,4);
z = vox2ras(3,4);

d = sqrt(sum(vox2ras(:,1:3).^2));
Mdc = vox2ras(1:3,1:3) ./ repmat(d,[3 1]);
if(det(Mdc) == 0.0)
  fprintf('ERROR: vox2ras determinant is 0\n');
  return;
end

r11 = Mdc(1,1);
r21 = Mdc(2,1);
r31 = Mdc(3,1);
r12 = Mdc(1,2);
r22 = Mdc(2,2);
r32 = Mdc(3,2);
r13 = Mdc(1,3);
r23 = Mdc(2,3);
r33 = Mdc(3,3);

if(det(Mdc) > 0.0) qfac = 1.0;
else
  r13 = -r13;
  r23 = -r23;
  r33 = -r33;
  qfac = -1.0;
end

%  /* following mat44_to_quatern() */

a = r11 + r22 + r33 + 1.0;
if(a > 0.5)
  a = 0.5 * sqrt(a);
  b = 0.25 * (r32-r23) / a;
  c = 0.25 * (r13-r31) / a;
  d = 0.25 * (r21-r12) / a;
else
  xd = 1.0 + r11 - (r22+r33);
  yd = 1.0 + r22 - (r11+r33);
  zd = 1.0 + r33 - (r11+r22);
  if(xd > 1.0)
    b = 0.5 * sqrt(xd);
    c = 0.25 * (r12+r21) / b;
    d = 0.25 * (r13+r31) / b;
    a = 0.25 * (r32-r23) / b;
  elseif( yd > 1.0 )
    c = 0.5 * sqrt(yd);
    b = 0.25 * (r12+r21) / c;
    d = 0.25 * (r23+r32) / c;
    a = 0.25 * (r13-r31) / c;
  else
    d = 0.5 * sqrt(zd);
    b = 0.25 * (r13+r31) / d;
    c = 0.25 * (r23+r32) / d;
    a = 0.25 * (r21-r12) / d;
  end
  if(a < 0.0)
    a = -a;
    b = -b;
    c = -c;
    d = -d;
  end
end


return
