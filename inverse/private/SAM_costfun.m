function inv_pseudo_Z = SAM_costfun(angle, position, tanu, tanv, lf, covariance, inv_covariance, noise_covariance)
        
% costfunction for non-linear beamformer. Use this cost-function to
% find the optimum orientation (in the tangential plane formed by
% tanu and tanv) of the targetvoxel maximizes the pseudo_Z (i.e.
% minimises the inverse of pseudo_Z)
%
% positions in mm in CTF co-ordinate system
%
% AH, 05april 2005: if origin = [], then the localspheres headmodel
% will be used for the forward calculations. The localspheres origins
% should be given in forward_resource (in mm in CTF co-ordinates)

MDip = settang(angle, tanu, tanv);
MagDip = sqrt(dot(MDip,MDip));
UnitMDip = MDip/MagDip;

surface.positions = position;
surface.vert_normals = UnitMDip;

gain = lf * UnitMDip';

% if ~isempty(origin)
%     % single sphere forward calculation
%     [junk, gain] = ComputeLeads(surface, forward_resource, origin, [], 0); 
% else
%     % multi sphere forward calculation
%     % disp(sprintf('\n\n USING LOCALSPHERES\n\n'))
%     gain = Compute_leads_Multi_fastM(forward_resource, surface.positions, surface.vert_normals, forward_resource.startsensor);
% end

trgain_invC = gain' * inv_covariance;
SAMweights =  trgain_invC / (trgain_invC * gain);

P = SAMweights * covariance * SAMweights';
N = SAMweights * noise_covariance * SAMweights'; 
inv_pseudo_Z = sqrt(N/P);

% make sure the angle stays between 0 and 180 degrees
if (angle > pi) || (angle < 0) 
    inv_pseudo_Z = inv_pseudo_Z*10^4;
else
%    hold on
%    plot(angle*180/pi, inv_pseudo_Z,'b*')
%    drawnow    
end  

