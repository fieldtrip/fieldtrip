%  HBF_TM_PHI_LC makes a BEM transfer matrix for potential using LC approach
% 
%  [T,startinds,endinds]=HBF_TM_PHI_LC(D,ci,co,zerolevel)
% 
%    D:  D matrix, cell [N(meshes) x N(meshes)]
%    ci: conductivity inside each boundary surface, [N(meshes) x 1]
%    co: conductivity outside each boundary surface, [N(meshes) x 1]
%    zerolevel (optional):  index to the mesh, on which the mean potential 
%        (over vertices) is set to zero; use 0 to omit 
%    T: BEM transfer matrix, N(vertices in the whole model)^2
%    startinds, endinds: indices that point to rows/columns of T that
%        correspond to each BEM boundary surface;
%        inds{1}=startinds(1):endinds(1)
%  
%    If the model is embedded in infinite conductor (that is, co(end)>0),
%    zero of the potential is uniquely defined and 'zerolevel' argument is
%    not needed. In a finite model (co(end)=0), the zero of the potential
%    needs to be specified. The default choice is that the mean potential 
%    over the vertices of the outer boundary of the model is zero.
% 
%  v160302 (c) Matti Stenroos
%