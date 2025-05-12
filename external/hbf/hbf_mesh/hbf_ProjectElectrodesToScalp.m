function elecs=hbf_ProjectElectrodesToScalp(elecpos,bmeshes)
% HBF_PROJECTELECTRODESTOSCALP projects a set of electrodes to the outermost
% mesh, makes some interpolation coefficients, and builds the hbf struct "elecs".
%
% elecs = HBF_PROJECTELECTRODESTOSCALP(elecpos,bmeshes)
% 
% elecpos:  positions of electrodes (somewhat coregistered...)
% bmeshes:  hbf struct bmeshes
%
% elecs:    hbf struct for electrodes
% v200629 (c) Matti Stenroos

elecs.p_orig=elecpos;
[elecs.p_proj,elecs.loctype,elecs.locinfo,elecs.projdist]=hbf_ProjectPointsToMesh(bmeshes{end},elecpos);
elecs.NtoE=sparse(hbf_LinearInterpolationMatrix(bmeshes{end},elecs.p_proj,elecs.loctype,elecs.locinfo));


