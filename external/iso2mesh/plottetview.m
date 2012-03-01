function plottetview(session,method)
%
% plottetview(session,method)
%
% wrapper for tetview to plot the generated mesh
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input:
%	 session: a string to identify the output files for plotting, '' for
%	          the default session
%    method:  method can be 'cgalsurf' (default), 'simplify', 'cgalpoly'
%             'cgalmesh' and 'remesh'
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

exesuff=getexeext;
exesuff=fallbackexeext(exesuff,'tetview');

sid=getvarfrom({'caller','base'},'ISO2MESH_SESSION');
if(nargin>=1)
  if(~ischar(session)) error('session name must be a string');  end
  if(~isempty(session))
     evalin('base',['ISO2MESH_SESSION=''' session ''';']);
  end
end

if(nargin<2) % method=cgalsurf by default
  cmd=sprintf('"%s%s" "%s"',mcpath('tetview'),exesuff,mwpath('post_vmesh.1'));
elseif(strcmp(method,'cgalmesh'))
  cmd=sprintf('"%s%s" "%s"',mcpath('tetview'),exesuff,mwpath('post_cgalmesh.mesh'));
elseif(strcmp(method,'cgalsurf'))
  cmd=sprintf('"%s%s" "%s"',mcpath('tetview'),exesuff,mwpath('post_extract.off'));
elseif(strcmp(method,'cgalpoly'))
  cmd=sprintf('"%s%s" "%s"',mcpath('tetview'),exesuff,mwpath('post_cgalpoly.mesh'));
elseif(strcmp(method,'remesh'))
  cmd=sprintf('"%s%s" "%s"',mcpath('tetview'),exesuff,mwpath('post_remesh.off'));
else
  if(~isempty(sid))
    evalin('base',['ISO2MESH_SESSION=''' sid ''';']);
  end
  error('unknown method');
end

if(~isempty(sid)) 
  evalin('base',['ISO2MESH_SESSION=''' sid ''';']); 
end

system(cmd);

