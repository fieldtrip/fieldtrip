function [mri, hdr] = read_ctf_mri(filename);

% READ_CTF_MRI reads header and imnage data from CTF format MRI file
%
% [mri, hdr] = read_ctf_mri(filename)
%
% See also READ_CTF_MEG4, READ_CTF_RES4

% Copyright (C) 2003 Robert Oostenveld
%
% $Log: read_ctf_mri.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.4  2008/09/30 07:47:04  roboos
% replaced all occurences of setstr() with char(), because setstr is deprecated by Matlab
%
% Revision 1.3  2005/08/26 13:49:03  roboos
% changed warp3d into warp_apply
%
% Revision 1.2  2003/07/23 15:02:27  roberto
% added check on valid input for read_ctf_meg4, other changes unknown
%
% Revision 1.1  2003/06/10 08:14:54  roberto
% new implementation
%

fid = fopen(filename,'rb', 'ieee-be');

if fid<=0
  error(sprintf('could not open MRI file: %s\n', filename));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ THE IMAGE HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off
% general header information
hdr.identifierString = char(fread(fid,32,'char'))'; % CTF_MRI_FORMAT VER 2.2
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
fread(fid,2,'char'); % padding to 4 byte boundary
hdr.HeadModel.defaultSphereX = fread(fid,1,'float'); % sphere origin x coordinate(in mm)
hdr.HeadModel.defaultSphereY = fread(fid,1,'float'); % sphere origin y coordinate(in mm)
hdr.HeadModel.defaultSphereZ = fread(fid,1,'float'); % sphere origin z coordinate(in mm)
hdr.HeadModel.defaultSphereRadius = fread(fid,1,'float'); % default sphere radius(in mm)

% Image_Info specific header items
hdr.Image.modality = fread(fid,1,'int16'); % 0 = MRI, 1 = CT, 2 = PET, 3 = SPECT, 4 = OTHER
hdr.Image.manufacturerName = char(fread(fid,64,'char'))';
hdr.Image.instituteName = char(fread(fid,64,'char'))';
hdr.Image.patientID = char(fread(fid,32,'char'))';
hdr.Image.dateAndTime = char(fread(fid,32,'char'))';
hdr.Image.scanType = char(fread(fid,32,'char'))';
hdr.Image.contrastAgent = char(fread(fid,32,'char'))';
hdr.Image.imagedNucleus = char(fread(fid,32,'char'))';
fread(fid,2,'char'); % padding to 4 byte boundary
hdr.Image.Frequency = fread(fid,1,'float');
hdr.Image.FieldStrength = fread(fid,1,'float');
hdr.Image.EchoTime = fread(fid,1,'float');
hdr.Image.RepetitionTime = fread(fid,1,'float');
hdr.Image.InversionTime = fread(fid,1,'float');
hdr.Image.FlipAngle = fread(fid,1,'float');
hdr.Image.NoExcitations = fread(fid,1,'int16');
hdr.Image.NoAcquisitions = fread(fid,1,'int16');
hdr.Image.commentString = char(fread(fid,256,'char'))';
hdr.Image.forFutureUse = char(fread(fid,64,'char'))';

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
transformMatrix = fread(fid,[4 4],'float')'; % transformation matrix head->MRI[column][row]
fread(fid,204,'char'); % unused, padding to 1028 bytes
warning on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ THE IMAGE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if hdr.dataSize==1
  mri = uint8(fread(fid, 256*256*256, 'uint8'));
elseif hdr.dataSize==2
  mri = uint16(fread(fid, 256*256*256, 'uint16'));
else
  error('unknown datasize in CTF mri file');
end
mri = reshape(mri, [256 256 256]);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO POST-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reorient the image data to obtain corresponding image data and transformation matrix
mri = permute(mri, [3 1 2]);	% this was determined by trial and error

% reorient the image data and the transformation matrix along the left-right direction
% remember that the fiducials in voxel coordinates also have to be flipped (see down)
mri = flipdim(mri, 1);
flip = [-1 0 0 256
         0 1 0 0
         0 0 1 0
         0 0 0 1    ];
transformMatrix = flip*transformMatrix;

% re-compute the homogeneous transformation matrices (apply voxel scaling)
scale = eye(4);
scale(1,1) = hdr.mmPerPixel_sagittal;
scale(2,2) = hdr.mmPerPixel_coronal;
scale(3,3) = hdr.mmPerPixel_axial;
hdr.transformHead2MRI = transformMatrix*inv(scale);
hdr.transformMRI2Head = scale*inv(transformMatrix);

% determint location of fiducials in MRI voxel coordinates
% flip the fiducials in voxel coordinates to correspond to the previous flip along left-right
hdr.fiducial.mri.nas = [256 - hdr.HeadModel.Nasion_Sag hdr.HeadModel.Nasion_Cor hdr.HeadModel.Nasion_Axi];
hdr.fiducial.mri.lpa = [256 - hdr.HeadModel.LeftEar_Sag hdr.HeadModel.LeftEar_Cor hdr.HeadModel.LeftEar_Axi];
hdr.fiducial.mri.rpa = [256 - hdr.HeadModel.RightEar_Sag hdr.HeadModel.RightEar_Cor hdr.HeadModel.RightEar_Axi];

% compute location of fiducials in MRI and HEAD coordinates
hdr.fiducial.head.nas = warp_apply(hdr.transformMRI2Head, hdr.fiducial.mri.nas, 'homogenous');
hdr.fiducial.head.lpa = warp_apply(hdr.transformMRI2Head, hdr.fiducial.mri.lpa, 'homogenous');
hdr.fiducial.head.rpa = warp_apply(hdr.transformMRI2Head, hdr.fiducial.mri.rpa, 'homogenous');

