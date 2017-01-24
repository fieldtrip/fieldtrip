function [output] = volumesmooth(input, fwhm, str)

% VOLUMESMOOTH is a helper function for segmentations
%
% See also VOLUMETHRESHOLD, VOLUMEFILLHOLES

if nargin==3
  fprintf('smoothing %s with a %d-voxel FWHM kernel\n', str, fwhm);
end

% ensure that SPM is available (any version)
ft_hastoolbox('spm12', 3) || ft_hastoolbox('spm8', 3) || ft_hastoolbox('spm2', 3);

% don't touch the input, make a deep copy
output = input+0;

% the mex files underneath the spm function will change the input variable
spm_smooth(output, output, fwhm);
