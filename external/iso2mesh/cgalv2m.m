function [node,elem,face]=cgalv2m(vol,opt,maxvol)
%
% [node,elem,face]=cgalv2m(vol,opt,maxvol)
%
% wrapper for CGAL 3D mesher (CGAL 3.5 or up)
% convert a binary (or multi-valued) volume to tetrahedral mesh
%
% http://www.cgal.org/Manual/3.5/doc_html/cgal_manual/Mesh_3/Chapter_main.html
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input:
%	 vol: a volumetric binary image
%	 ix,iy,iz: subvolume selection indices in x,y,z directions
%	 opt: parameters for CGAL mesher, if opt is a structure, then
%	     opt.radbound: defines the maximum surface element size
%	     opt.angbound: defines the miminum angle of a surface triangle
%	     opt.distbound: defines the maximum distance between the 
%		 center of the surface bounding circle and center of the 
%		 element bounding sphere
%	     opt.reratio:  maximum radius-edge ratio
%	     if opt is a scalar, it only specifies radbound.
%	 maxvol: target maximum tetrahedral elem volume
%
% output:
%	 node: output, node coordinates of the tetrahedral mesh
%	 elem: output, element list of the tetrahedral mesh, the last 
%	      column is the region id
%	 face: output, mesh surface element list of the tetrahedral mesh
%	      the last column denotes the boundary ID
%	      note: each triangle will appear twice in the face list with each
%		    one attaches to each side of the interface. one can remove
%		    the redundant triangles by unique(face(:,1:3),'rows')
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

fprintf(1,'creating surface and tetrahedral mesh from a multi-domain volume ...\n');

dtype=class(vol);
if(~(islogical(vol) | strcmp(dtype,'uint8')))
	error('cgalmesher can only handle uint8 volumes, you have to convert your image to unit8 first.');
end

if(~any(vol))
        error('no labeled regions found in the input volume.');
end

exesuff=getexeext;
exesuff=fallbackexeext(exesuff,'cgalmesh');

ang=30;
ssize=6;
approx=0.5;
reratio=3;

if(~isstruct(opt))
	ssize=opt;
end

if(isstruct(opt) & length(opt)==1)  % does not support settings for multiple labels
	if(isfield(opt,'radbound'))   ssize=opt.radbound; end
	if(isfield(opt,'angbound'))   ang=opt.angbound; end
	if(isfield(opt,'distbound')) approx=opt.distbound; end
	if(isfield(opt,'reratio'))    reratio=opt.reratio; end
end

saveinr(vol,mwpath('pre_cgalmesh.inr'));
deletemeshfile(mwpath('post_cgalmesh.mesh'));

randseed=hex2dec('623F9A9E'); % "U+623F U+9A9E"

if(~isempty(getvarfrom('base','ISO2MESH_RANDSEED')))
        randseed=getvarfrom('base','ISO2MESH_RANDSEED');
end

cmd=sprintf('"%s%s" "%s" "%s" %f %f %f %f %f %d',mcpath('cgalmesh'),exesuff,...
    mwpath('pre_cgalmesh.inr'),mwpath('post_cgalmesh.mesh'),ang,ssize,...
    approx,reratio,maxvol,randseed);
system(cmd);
if(~exist(mwpath('post_cgalmesh.mesh'),'file'))
    error(['output file was not found, failure was encountered when running command: \n',cmd]);
end
[node,elem,face]=readmedit(mwpath('post_cgalmesh.mesh'));

% if a transformation matrix/offset vector supplied, apply them
if (isstruct(opt) & length(opt)==1)
    if(isfield(opt,'A') & isfield(opt,'B'))
        node(:,1:3)=(opt.A*node(:,1:3)'+repmat(opt.B(:),1,size(node,1)))';
    end
end

fprintf(1,'node number:\t%d\ntriangles:\t%d\ntetrahedra:\t%d\nregions:\t%d\n',...
    size(node,1),size(face,1),size(elem,1),length(unique(elem(:,end))));
fprintf(1,'surface and volume meshes complete\n');

