function [node,elem]=vol2restrictedtri(vol,thres,cent,brad,ang,radbound,distbound,maxnode)
%
% [node,elem]=vol2restrictedtri(vol,thres,cent,brad,ang,radbound,distbound,maxnode)
%
% surface mesh extraction using CGAL mesher
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2009/01/06
%
% input:
%       vol: a 3D volumetric image
%       thres: a scalar as the threshold of of the extraction
%       cent: a 3d position (x,y,z) which locates inside the resulting
%             mesh, this is automatically computed from vol2surf
%       brad: maximum bounding sphere squared of the resulting mesh
%       ang: minimum angular constrains of the resulting tranglar elements
%            (in degrees)
%       radbound: maximum triangle delaunay circle radius
%       distbound: maximum delaunay sphere distances
%       maxnode: maximum number of surface nodes (even radbound is not reached)
% output:
%       node: the list of 3d nodes in the resulting surface (x,y,z)
%       elem: the element list of the resulting mesh (3 columns of integers)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(radbound<1)
    warning(['You are meshing the surface with sub-pixel size. If this ' ...
             'is not your your intent, please check if you set ' ...
             '"opt.radbound" correctly for the default meshing method.']);
end

exesuff=getexeext;
exesuff=fallbackexeext(exesuff,'cgalsurf');

saveinr(vol,mwpath('pre_extract.inr'));
deletemeshfile(mwpath('post_extract.off'));

randseed=hex2dec('623F9A9E'); % "U+623F U+9A9E"

if(~isempty(getvarfrom({'caller','base'},'ISO2MESH_RANDSEED')))
	randseed=getvarfrom({'caller','base'},'ISO2MESH_RANDSEED');
end

initnum=50;
if(~isempty(getvarfrom({'caller','base'},'ISO2MESH_INITSIZE')))
        initnum=getvarfrom({'caller','base'},'ISO2MESH_INITSIZE');
end

system([' "' mcpath('cgalsurf') exesuff '" "' mwpath('pre_extract.inr') ...
    '" ' sprintf('%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %d ',thres,cent,brad,ang,radbound,distbound,maxnode) ...
    ' "' mwpath('post_extract.off') '" ' sprintf('%.0f %d',randseed,initnum)]);
[node,elem]=readoff(mwpath('post_extract.off'));

% assuming the origin [0 0 0] is located at the lower-bottom corner of the image
node=node+0.5;
