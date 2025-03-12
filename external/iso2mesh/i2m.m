function varargout = i2m(varargin)
%
% newworkspace=i2m;
%  or
% newworkspace=i2m(workspace)
%
% Shortcut for img2mesh, a GUI for iso2mesh
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input/output: please see details in the help for img2mesh
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

[varargout{1:nargout}] = img2mesh(varargin{:});
