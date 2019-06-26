function [lf] = leadfield_simbio(pos, vol)

% leadfield_simbio leadfields for a set of dipoles
%
% [lf] = leadfield_simbio(pos, vol);
%
% with input arguments
%   pos     a matrix of dipole positions
%           there can be 'deep electrodes' too!
%   vol     contains a FE volume conductor (output of ft_prepare_vol_sens)
%
% the output lf is the leadfield matrix of dimensions m (rows) x n*3 (cols)

% copyright (c) 2012, Johannes Vorwerk

try   
    lf = zeros(3*size(pos,1),size(vol.transfer,1));
    dir = diag([1,1,1]);
    for i=1:size(pos,1)
        locpos = repmat(pos(i,:),3,1);
        rhs = sb_rhs_venant(locpos,dir,vol);
        lf((3*(i-1)+1):(3*(i-1)+3),:) = (vol.transfer * rhs)';
    end    
    lf = lf';
catch
  ft_warning('an error occurred while running simbio');
  rethrow(lasterror)
end
