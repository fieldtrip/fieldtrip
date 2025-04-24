%  HBF_LFM_BP builds lead field matrix for primary magnetic field. 
%  conductor
% 
%  LFM=HBF_LFM_BP(coils,spos,sdir)
%  LFM=HBF_LFM_BP(coils,spos)
% 
%    coils:  coil description, hbf struct
%    spos:   source positions, [M x 3]
%    sdir:   source orientations (unit-length) or moments, [M x 3]
% 
%    LFM:   lead field matrix, [Number of coils (field points) x M]
%            [t_1 ... t_M]
%            or, if 'sdir' omitted, [Number of coils x 3M]
%            [t_1x t_1y t1_z ... t_Mx t_My t_Mz]
% 
%  You can also compute primary magnetic field due to any set of directed
%  dipoles by giving the dipole moments (with amplitude) in the 'sdir'
%  argument. If you omit 'sdir', LFM will be computed using xyz unit dipole
%  triplets
% 
%  "Primary magnetic field" B_primary is the magnetic field that arises from
%  the source current only, without contribution from the volume conductor.
%  It is the same as the magnetic field due to the source in an infinite,
%  homogeneous volume conductor. The total magnetic field as computed by
%  HBF_LFM_B is
%  B_total = B_primary + B_volume,
%  where the latter term shows the contribution of volume currents driven by
% the electric field.
% 
%  v201022 (c) Matti Stenroos
%