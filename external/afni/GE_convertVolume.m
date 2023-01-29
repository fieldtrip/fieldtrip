function status = GE_convertVolume(inDir,runnum,outName)
%
% status = GE_convertVolume(inDir,runnum,outName)
%
%GE_convertVolume
%
% Version 3.1
%
% Converts a series of GE slices into Analyze format for SPM.
% inDir is the name of the directory containing the first file of
% the series, e.g. 003
% runnum is the run number of the set you want to convert
% outName is used for naming the Analyze files
% status = 1 if there is an error.
%
% Modified to use spm functions
%
% Souheil J. Inati
% Dartmouth College
% September 2001
% souheil.inati@dartmouth.edu
%

status = 1;

if (nargin < 3)
        error('Not enough input arguments.')
        return
end

% Create the name of the first file in inDir
firstfile = fullfile(inDir,'I.001');

% Read the Header from the first file
[su_hdr,ex_hdr,se_hdr,im_hdr,pix_hdr,im_offset] = GE_readHeader(firstfile);

% The image Dimensions
Nx =  im_hdr.imatrix_X; % X Voxels
Ny =  im_hdr.imatrix_Y; % Y Voxels
Nz =  im_hdr.slquant;   % Z Voxels
volSize = [Nx Ny Nz];

% Read in the volume
[imageVol, lastfile] = GE_readVolume(inDir, runnum, volSize, pix_hdr.img_depth, im_offset);

% Is the first image the first or the last?
if (se_hdr.end_loc - se_hdr.start_loc) > 0
  scandir = 1;
else
  scandir = -1;
end

% Compute the M matrix
M = GE_createSPMmat(im_hdr,scandir);

% Create the SPM volume
V.fname = outName;
V.dim = [volSize spm_type('int16')]; % short ints
V.mat = M;
V.pinfo = [1 0 0]';
V = spm_create_image(V);

% Write out the SPM Volume
% Don't use spm_write_vol because that sets a scale factor.
for i = 1:Nz
  V = spm_write_plane(V,squeeze(imageVol(:,:,i)),i);
end

% Done with vol
fprintf('Wrote %s.img.\n',outName);

% Done
fprintf('\nConversion Finished. \n');

status = 0;

return


