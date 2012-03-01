function openface=volface(t)
%
% openface=volface(t)
%
% find the surface patches of a volume
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2009/10/13
%
% input:
%      t: input, volumetric element list, dimension (ne,4)
%
% output:
%      openface: list of faces of the specified volume
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

openface=surfedge(t);
