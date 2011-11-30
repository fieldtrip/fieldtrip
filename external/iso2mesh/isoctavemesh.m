function isoctave=isoctavemesh
%
% isoctave=isoctavemesh
%
% determine whether the code is running in octave
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% output:
%   isoctave: 1 if in octave, otherwise 0
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%
isoctave=(exist('OCTAVE_VERSION')~=0);
