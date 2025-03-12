function varargout = volface(t)
%
% [openface,elemid]=volface(t)
%
% find the surface patches of a volume
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2009/10/13
%
% input:
%      t: input, volumetric element list, dimension (ne,4)
%
% output:
%      openface: list of faces of the specified volume
%      elemid (optional): the corresponding index of the
%                tetrahedron of an open-edge or triangle,
%                elemid has the same length as openedge.
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

[varargout{1:nargout}] = surfedge(t);
