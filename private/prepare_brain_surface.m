function [pnt] = prepare_brain_surface(cfg, grad, vol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that makes a brain surface from the headhape and places a lot
% of dipoles on it. The dipoles subsequently can be used for a simplified 
% distributed source model.
%
% This function uses the following fields from the configuration
%   cfg.headshape
%   cfg.spheremesh
%   cfg.inwardshift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2004, Robet Oostenveld
%
% $Log: prepare_brain_surface.m,v $
% Revision 1.6  2007/08/06 09:20:14  roboos
% added support for bti_hs
%
% Revision 1.5  2005/11/16 09:10:47  roboos
% only change in some whitespace
%
% Revision 1.4  2005/11/04 22:16:59  geekra
% Changed triangulation method when constructing surface from multiple sphere
% model. Now all points are used to create the surface, preventing problems
% with the function normals.
%
% Revision 1.3  2005/11/02 10:09:21  geekra
% Removed extra bracket in line 73 that resulted in a syntax error.
%
% Revision 1.2  2005/11/01 09:52:00  roboos
% Changed the construction of the dipole mesh for a localspheres
% headmodel, which was hardcoded for the 151 channel CTF system. Now
% it uses an updated version of the head_surf function. Furthermore,
% implemented support for cfg.spheremesh values other than the number
% of vertices of the refined icosahedrons, and made the reducepatch
% of the headshape dependent on the cfg.spheremesh.
%
% Revision 1.1  2004/06/28 08:59:38  roboos
% moved files from fieldtrip to fieldtrip/private
%
% Revision 1.3  2004/05/17 07:22:16  roberto
% fixed bug in triangle orientation for multi-sphere
%
% Revision 1.2  2004/04/26 12:36:49  roberto
% fixed bug in case of multisphere head surface
%
% Revision 1.1  2004/04/08 15:51:20  roberto
% initial submissiion into cvs
%

% create a head surface, which by default is assumed to correspond to the skin
if strcmp(cfg.headshape, 'headmodel')
  if length(vol.r)==1
    % single sphere model, use a refined icosahedron 
    fprintf('constructing brain surface from single sphere model\n');
    if cfg.spheremesh==2562
      [pnt, tri] = icosahedron2562;
    elseif cfg.spheremesh==642
      [pnt, tri] = icosahedron642;
    elseif cfg.spheremesh==162
      [pnt, tri] = icosahedron162;
    elseif cfg.spheremesh==42
      [pnt, tri] = icosahedron42;
    elseif cfg.spheremesh==12
      [pnt, tri] = icosahedron;
    else
      [pnt, tri] = ksphere(cfg.spheremesh);
    end
    % scale the sourcemodel sphere to the size of the headmodel sphere
    pnt = pnt * vol.r;
    pnt(:,1) = pnt(:,1) + vol.o(1);
    pnt(:,2) = pnt(:,2) + vol.o(2);
    pnt(:,3) = pnt(:,3) + vol.o(3);
  else
    % multiple sphere model, use the points on the skin surface
    fprintf('constructing brain surface from multiple sphere model\n');
    [pnt, tri] = head_surf(vol, grad, 0);
    prj = elproj(pnt);
    tri = delaunay(prj(:,1),prj(:,2));
  end
else
  fprintf('constructing brain surface from headshape file\n');
  % read the headshape from file
  if filetype(cfg.headshape, 'ctf_shape')
    shape = read_ctf_shape(cfg.headshape);
    pnt = shape.pnt;
  elseif filetype(cfg.headshape, '4d_hs')
    pnt = read_bti_hs(cfg.headshape);
  end
  prj = elproj(pnt);
  tri = delaunay(prj(:,1), prj(:,2));
  % the number of triangles is approximately twice the number of vertices
  [tri, pnt] = reducepatch(tri, pnt, 2*cfg.spheremesh);
end

% correct the orientation of the triangles
sel=find(solid_angle(pnt,tri)<0);
tri(sel,:) = fliplr(tri(sel,:));

% shift the head surface inward with a certain amount
ori = normals(pnt, tri, 'vertex');
pnt = pnt + cfg.inwardshift * ori;
