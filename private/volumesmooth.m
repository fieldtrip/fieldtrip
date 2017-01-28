function [output] = volumesmooth(input, fwhm, str)

% VOLUMESMOOTH is a helper function for segmentations
%
% See also VOLUMETHRESHOLD, VOLUMEFILLHOLES

% ensure that SPM is available, needed for spm_smooth
hasspm = ft_hastoolbox('spm8up', 3) || ft_hastoolbox('spm2', 1);

if nargin==3
  fprintf('smoothing %s with a %d-voxel FWHM kernel\n', str, fwhm);
end

% don't touch the input, make a deep copy
output = input+0;

% the mex files underneath the spm function will change the input variable
spm_smooth(output, output, fwhm);
