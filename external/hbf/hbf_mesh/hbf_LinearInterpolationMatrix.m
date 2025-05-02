%  HBF_LINEARINTERPOLATIONMATRIX makes a matrix for interpolating a function
%  in electrode positions
% 
%  function [ntoe_full,nodeinds,ntoe_dense]=LinearInterpolationMatrix(mesh,e_proj,e_loctype,e_locinfo)
% 
%  mesh:       boundary mesh, on which the electrodes are projected
%  e_proj, e_loctype, e_loctype: electrode projection info from hbf struct
%              elecs, as obtained with HBF_PROJECTELECTRODESTOSCALP
% 
%  ntoe_full:  [N_elec x mesh.nop] matrix for interpolating a function
%              from full mesh onto projected electrode locations 
%  nodeinds:   [N x 1]; nodes that are actually needed for the electrodes
%  ntoe_dense: [N_elec x N]; the same as ntoe_full but zero-columns removed
% 
%  use of ntoe_dense:
%  T_temp=Tmesh(nodeinds,:);
%  T_elec=ntoe_dense*T_temp;
% 
%  v200924 (c) Matti Stenroos
%