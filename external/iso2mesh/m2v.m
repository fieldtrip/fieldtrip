function varargout=m2v(varargin)
%
% vol=m2v(node,face,Nxyz)
%  or
% vol=m2v(node,face,xi,yi,zi)
%
% shortcut for mesh2vol, rasterizing a teterahedral mesh to a volume using graphics
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input/output: please see details in the help for mesh2vol
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

[varargout{1:nargout}]=mesh2vol(varargin{:});
