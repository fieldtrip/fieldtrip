function faces = meshface(elem, varargin)
%
% faces=meshface(elem,opt)
%
% return all faces in a surface or volumetric mesh
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2011/02/26
%
% input:
%    elem:  element table of a mesh (support N-d space element)
%    opt: optional input, giving the additional options. If opt
%         is a struct, it can have the following field:
%       opt.nodeorder: if 1, assuming the elem node indices is in CCW
%                      orientation; 0 use nchoosek() output to order faces
%         you can replace opt by a series of ('param', value) pairs.
%
% output:
%    face:  face list; each row is an face, specified by the starting and
%           ending node indices, the total face number is
%           size(elem,1) x nchoosek(size(elem,2),2). All faces are ordered
%           by looping through each element first.
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

dim = size(elem);
faceid = nchoosek(1:dim(2), 3);
len = size(faceid, 1);
faces = zeros(dim(1) * len, 3);
for i = 0:len - 1
    faces((i * dim(1) + 1):((i + 1) * dim(1)), :) = [elem(:, faceid(i + 1, 1)) elem(:, faceid(i + 1, 2)) elem(:, faceid(i + 1, 3))];
end
