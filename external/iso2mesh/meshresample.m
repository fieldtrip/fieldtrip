function [node,elem]=meshresample(v,f,keepratio)
%
% [node,elem]=meshresample(v,f,keepratio)
%
% resample mesh using CGAL mesh simplification utility
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2007/11/12
%
% input:
%    v: list of nodes
%    f: list of surface elements (each row for each triangle)
%    keepratio: decimation rate, a number less than 1, as the percentage
%               of the elements after the sampling
%
% output:
%    node: the node coordinates of the sampled surface mesh
%    elem: the element list of the sampled surface mesh
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

[node,elem]=domeshsimplify(v,f,keepratio);

if(length(node)==0)
    warning(['Your input mesh contains topological defects, and the ',...
           'mesh resampling utility aborted during processing. Now iso2mesh ',...
           'is trying to repair your mesh with meshcheckrepair. ',...
           'You can also call this manually before passing your mesh to meshresample.'] );
    [vnew,fnew]=meshcheckrepair(v,f);
    [node,elem]=domeshsimplify(vnew,fnew,keepratio);
end
[node,I,J]=unique(node,'rows');
elem=J(elem);
saveoff(node,elem,mwpath('post_remesh.off'));

end

% function to perform the actual resampling
function [node,elem]=domeshsimplify(v,f,keepratio)
  exesuff=getexeext;
  exesuff=fallbackexeext(exesuff,'cgalsimp2');

  saveoff(v,f,mwpath('pre_remesh.off'));
  deletemeshfile(mwpath('post_remesh.off'));
  system([' "' mcpath('cgalsimp2') exesuff '" "' mwpath('pre_remesh.off') '" ' num2str(keepratio) ' "' mwpath('post_remesh.off') '"']);
  [node,elem]=readoff(mwpath('post_remesh.off'));
end
