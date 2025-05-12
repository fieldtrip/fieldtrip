%  HBF_TM_PHI_LC_ISA2 makes a BEM transfer matrix for potential using LC 
%    BEM and isolated source approach
% 
%  [T,startinds,endinds]=HBF_TM_PHI_LC_ISA2(D,ci,co,ISAsurf)
%    D:  D matrix, cell [N(meshes) x N(meshes)]
%    ci: conductivity inside each boundary surface, [N(meshes) x 1]
%    co: conductivity outside each boundary surface, [N(meshes) x 1]
%    ISAsurf:  index to the mesh, on which the isolation is performed. 
% 
%    T: BEM transfer matrix, N(vertices in the whole model)^2
%    startinds, endinds: indices that point to rows/columns of T that
%        correspond to each BEM boundary surface;
%        inds{1}=startinds(1):endinds(1)
%  
%    In EEG/MEG application, the isolation is typically done on the inner skull
%    boundary (mesh index 1 in three-shell model, 2 or 3 in a four-compartment
%    model). All sources must be inside the isolation boundary.
% 
%    The mean potential over the nodes of the isolation surface is set to zero.
%    Stenroos & Sarvas, PMB 2012
%  v160303 (c) Matti Stenroos
%