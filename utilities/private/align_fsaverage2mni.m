function [mri] = align_fsaverage2mni(mri)

% ALIGN_FSAVERAGE2MNI performs an affine alignment of the anatomical volume from
% FSAVERAGE towards MNI coordinates. Only the homogeneous transformation matrix is
% modified and the coordsys-field is updated.
%
% Use as
%   mri = align_fsaverage2mni(mri)
% where the first input argument is a FieldTrip MRI-structure.
%
% with fsaverage we mean MNI305
% with mni       we mean MNI152, i.e. the template used in SPM
%
% See http://freesurfer.net/fswiki/CoordinateSystems

fsaverage2mni = [
   0.9975   -0.0073    0.0176   -0.0429
   0.0146    1.0009   -0.0024    1.5496
  -0.0130   -0.0093    0.9971    1.1840
  ];

assert(strcmp(mri.coordsys, 'fsaverage'), 'incorrect input coordinate system ''%s''', mri.coordsys);
mri.transform = fsaverage2mni * mri.transform;
mri.coordsys = 'mni';
