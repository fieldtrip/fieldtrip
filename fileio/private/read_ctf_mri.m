function [mri, hdr] = read_ctf_mri(filename)


% READ_CTF_MRI reads header and image data from a CTF version 2.2 MRI file
%
% Use as
%   [mri, hdr] = read_ctf_mri(filename)
%
% See also READ_CTF_MRI4

% Copyright (C) 2003-2010 Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% Some versions require specifying latin1 (ISO-8859-1) character encoding.
fid = fopen(filename, 'rb', 'ieee-be', 'ISO-8859-1');

if fid<=0
  error(sprintf('could not open MRI file: %s\n', filename));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ THE IMAGE HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ws = warning('off');

% general header information
hdr.identifierString = fread(fid,[1 32],'uint8=>char'); % CTF_MRI_FORMAT VER 2.2
hdr.imageSize = fread(fid,1,'int16'); % always = 256
hdr.dataSize = fread(fid,1,'int16'); % 1 or 2(bytes)
hdr.clippingRange = fread(fid,1,'int16'); % max.integer value of data
hdr.imageOrientation = fread(fid,1,'int16'); % eg., 0 = left on left, 1 = left on right
hdr.mmPerPixel_sagittal = fread(fid,1,'float'); % voxel dimensions in mm
hdr.mmPerPixel_coronal = fread(fid,1,'float'); % voxel dimensions in mm
hdr.mmPerPixel_axial = fread(fid,1,'float'); % voxel dimensions in mm

% HeadModel_Info specific header items
hdr.HeadModel.Nasion_Sag = fread(fid,1,'int16'); % fid.point coordinate(in voxels) for nasion - sagittal
hdr.HeadModel.Nasion_Cor = fread(fid,1,'int16'); % nasion - coronal
hdr.HeadModel.Nasion_Axi = fread(fid,1,'int16'); % nasion - axial
hdr.HeadModel.LeftEar_Sag = fread(fid,1,'int16'); % left ear - sagittal
hdr.HeadModel.LeftEar_Cor = fread(fid,1,'int16'); % left ear - coronal
hdr.HeadModel.LeftEar_Axi = fread(fid,1,'int16'); % left ear - axial
hdr.HeadModel.RightEar_Sag = fread(fid,1,'int16'); % right ear - sagittal
hdr.HeadModel.RightEar_Cor = fread(fid,1,'int16'); % right ear - coronal
hdr.HeadModel.RightEar_Axi = fread(fid,1,'int16'); % right ear - axial
fread(fid,2,'uint8'); % padding to 4 byte boundary
hdr.HeadModel.defaultSphereX = fread(fid,1,'float'); % sphere origin x coordinate(in mm)
hdr.HeadModel.defaultSphereY = fread(fid,1,'float'); % sphere origin y coordinate(in mm)
hdr.HeadModel.defaultSphereZ = fread(fid,1,'float'); % sphere origin z coordinate(in mm)
hdr.HeadModel.defaultSphereRadius = fread(fid,1,'float'); % default sphere radius(in mm)

% Image_Info specific header items
hdr.Image.modality = fread(fid,1,'int16'); % 0 = MRI, 1 = CT, 2 = PET, 3 = SPECT, 4 = OTHER
hdr.Image.manufacturerName = fread(fid,[1 64],'uint8=>char');
hdr.Image.instituteName = fread(fid,[1 64],'uint8=>char');
hdr.Image.patientID = fread(fid,[1 32],'uint8=>char');
hdr.Image.dateAndTime = fread(fid,[1 32],'uint8=>char');
hdr.Image.scanType = fread(fid,[1 32],'uint8=>char');
hdr.Image.contrastAgent = fread(fid,[1 32],'uint8=>char');
hdr.Image.imagedNucleus = fread(fid,[1 32],'uint8=>char');
fread(fid,2,'uint8'); % padding to 4 byte boundary
hdr.Image.Frequency = fread(fid,1,'float');
hdr.Image.FieldStrength = fread(fid,1,'float');
hdr.Image.EchoTime = fread(fid,1,'float');
hdr.Image.RepetitionTime = fread(fid,1,'float');
hdr.Image.InversionTime = fread(fid,1,'float');
hdr.Image.FlipAngle = fread(fid,1,'float');
hdr.Image.NoExcitations = fread(fid,1,'int16');
hdr.Image.NoAcquisitions = fread(fid,1,'int16');
hdr.Image.commentString = fread(fid,[1 256],'uint8=>char');
hdr.Image.forFutureUse = fread(fid,[1 64],'uint8=>char');

% continuation general header
hdr.headOrigin_sagittal = fread(fid,1,'float'); % voxel location of head origin
hdr.headOrigin_coronal = fread(fid,1,'float'); % voxel location of head origin
hdr.headOrigin_axial = fread(fid,1,'float'); % voxel location of head origin
% euler angles to align MR to head coordinate system(angles in degrees !)
hdr.rotate_coronal = fread(fid,1,'float'); % 1. rotate in coronal plane by this angle
hdr.rotate_sagittal = fread(fid,1,'float'); % 2. rotate in sagittal plane by this angle
hdr.rotate_axial = fread(fid,1,'float'); % 3. rotate in axial plane by this angle
hdr.orthogonalFlag = fread(fid,1,'int16'); % if set then image is orthogonal
hdr.interpolatedFlag = fread(fid,1,'int16'); % if set than image was interpolated
hdr.originalSliceThickness = fread(fid,1,'float'); % original spacing between slices before interpolation
hdr.transformMatrix = fread(fid,[4 4],'float')'; % transformation matrix head->MRI[column][row]

% Go to image data file position. 
% fread(fid,202,'uint8'); % unused, padding to 1028 bytes.
% The previous (commented) line should read 202 bytes to get to the end of
% the header (position 1028), but it seems some versions of Matlab (or
% perhaps only on some systems) doesn't read 2 bytes somewhere and end up
% in position 1026...  In any case, it caused an error with some files so
% we must explicitely seek to position 1028.
fseek(fid, 1028, 'bof');

% turn all warnings back on
warning(ws);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ THE IMAGE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if hdr.dataSize == 1
  precision = '*uint8';
elseif hdr.dataSize == 2
  if hdr.clippingRange < 2^15
    % I think this is usually the case, i.e. data is stored as signed 16
    % bit int, even though data is only positive.
    precision = '*int16';
  else
    precision = '*uint16';
  end
else
  error('unknown datasize (%d) in CTF mri file.', hdr.dataSize);
end
mri = fread(fid, hdr.imageSize.^3, precision);
mri = reshape(mri, [hdr.imageSize hdr.imageSize hdr.imageSize]);
fclose(fid);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO POST-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data is stored in PIR order (i.e. fastest changing direction goes from
% anterior to Posterior, then from superior to Inferior, finally from left
% to Right) assuming subject position was not too oblique and proper
% conversion to .mri format.  (Does this depend on hdr.imageOrientation?)
% On the other hand, the transformation matrix and fiducials are in RPI
% order so we must reorient the image data to match them.
mri = permute(mri, [3 1 2]);

transformMatrix = hdr.transformMatrix;

% Construct minimal transformation matrix if fiducials were not defined.  
% Bring data into same orientation (head coordinates are ALS) at least.
if all(transformMatrix == 0)
  transformMatrix(1, 2) = -1;
  transformMatrix(2, 1) = -1;
  transformMatrix(3, 3) = -1;
  transformMatrix(4, 4) = 1;
end

% determine location of fiducials in MRI voxel coordinates
% flip the fiducials in voxel coordinates to correspond to the previous flip along left-right
hdr.fiducial.mri.nas = [hdr.HeadModel.Nasion_Sag hdr.HeadModel.Nasion_Cor hdr.HeadModel.Nasion_Axi];
hdr.fiducial.mri.lpa = [hdr.HeadModel.LeftEar_Sag hdr.HeadModel.LeftEar_Cor hdr.HeadModel.LeftEar_Axi];
hdr.fiducial.mri.rpa = [hdr.HeadModel.RightEar_Sag hdr.HeadModel.RightEar_Cor hdr.HeadModel.RightEar_Axi];

% Reorient the image data, the transformation matrix and the fiducials
% along the left-right direction.
% This may have been done only for visualization?  It can probably be
% "turned off" without problem.
if false
  mri = flipdim(mri, 1);
  flip = [-1 0 0 hdr.imageSize+1
           0 1 0 0
           0 0 1 0
           0 0 0 1    ];
  transformMatrix = flip*transformMatrix;
  hdr.fiducial.mri.nas = [hdr.fiducial.mri.nas, 1] * flip(1:3, :)';
  hdr.fiducial.mri.lpa = [hdr.fiducial.mri.lpa, 1] * flip(1:3, :)';
  hdr.fiducial.mri.rpa = [hdr.fiducial.mri.rpa, 1] * flip(1:3, :)';
end

% re-compute the homogeneous transformation matrices (apply voxel scaling)
scale = eye(4);
scale(1,1) = hdr.mmPerPixel_sagittal;
scale(2,2) = hdr.mmPerPixel_coronal;
scale(3,3) = hdr.mmPerPixel_axial;
hdr.transformHead2MRI = transformMatrix*inv(scale);
hdr.transformMRI2Head = scale*inv(transformMatrix);

% compute location of fiducials in MRI and HEAD coordinates
hdr.fiducial.head.nas = ft_warp_apply(hdr.transformMRI2Head, hdr.fiducial.mri.nas, 'homogenous');
hdr.fiducial.head.lpa = ft_warp_apply(hdr.transformMRI2Head, hdr.fiducial.mri.lpa, 'homogenous');
hdr.fiducial.head.rpa = ft_warp_apply(hdr.transformMRI2Head, hdr.fiducial.mri.rpa, 'homogenous');


