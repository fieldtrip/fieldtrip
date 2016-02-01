function snorm=surfnorm(node,face,varargin)
%
% snorm=surfnorm(node,face)
%    or
% snorm=surfnorm(node,face,'Normalize',0)
%
% compute the normal vectors for a triangular surface
%
% author: Qianqian Fang, <fangq at nmr.mgh.harvard.edu>
%
% input:
%   node: a list of node coordinates (nn x 3)
%   face: a surface mesh triangle list (ne x 3)
%   opt: a list of optional parameters, currently surfnorm supports:
%        'Normalize': [1|0] if set to 1, the normal vectors will be 
%                           unitary (default)
%
% output:
%   snorm: output surface normal vector at each face
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

opt=varargin2struct(varargin{:});

snorm=surfplane(node,face);
snorm=snorm(:,1:3);

if(getoptkey('Normalize',1,opt))
    snorm=snorm./repmat(sqrt(sum(snorm.*snorm,2)),1,3);
end
