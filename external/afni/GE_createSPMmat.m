%%%%%%%%%%%%%%%%%%%%%%%%
%                      %
% Write SPM mat file   %
%                      %
%%%%%%%%%%%%%%%%%%%%%%%%
function M = GE_createSPMmat(im_hdr, scandir)
%
% Put the appropriate translations and rotations to the M matrix
% given the information in the image header and the direction of
% acquisition
%
% S. Inati
% Dartmouth College
% Apr. 2001
%

% The conversion from pixels to mm
Dims = diag( [im_hdr.pixsize_X, ...
	      im_hdr.pixsize_Y, ...
	      im_hdr.slthick + im_hdr.scanspacing ]);

% Compute the coordinate system in the image plane
tlhc = [ im_hdr.tlhc_R; im_hdr.tlhc_A; im_hdr.tlhc_S ]; % Top Left Hand Corner of Image
trhc = [ im_hdr.trhc_R; im_hdr.trhc_A; im_hdr.trhc_S ]; % Top Right Hand Corner of Image
brhc = [ im_hdr.brhc_R; im_hdr.brhc_A; im_hdr.brhc_S ]; % Bottom Right Hand Corner of Image

x = trhc - tlhc; x = x./sqrt(x'*x);   % xhat
y = trhc - brhc; y = y./sqrt(y'*y);   % yhat

% The normal to the plane
norm = [ im_hdr.norm_R; im_hdr.norm_A; im_hdr.norm_S ];
% The directional normal
z = scandir * norm;   % zhat

% Build the M matrix for SPM
% M takes a voxel from the image and gives it a coordinate in mm
% On the scanner the voxels start in the top left hand corner of the
% first image.  In SPM they start in the bottom left hand corner,
% so flip y and set the origin to tlhc.
% NB: The voxel (1,1,1) should have position tlhc
Rot = [x, -y, z];
M = eye(4);
M(1:3,1:3) = Rot  * Dims;
M(1:3,4) = tlhc - Rot*Dims*[1;1;1];

return