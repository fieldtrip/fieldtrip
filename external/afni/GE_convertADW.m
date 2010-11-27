function status = GE_convertADW(inDir,outStem,starttime,nvols)
%
% status = GE_convertADW(inDir,outStem,starttime,nvols)
%
%GE_convertADW
%
% Version 3.1
%
% Converts a series of GE slices acquired on the Advanced
% Development Workstation into Analyze format for SPM.
% inDir is the name of the first directory containing the
%      series, e.g. 003
% outStem is the stem used for naming the Analyze files
% starttime is the time point number to start
% nvols is the number of volumes (time points) to convert
%
% status = 1, error
% status = 0, all is well
%
% eg GE_convert('DATA/00023/003','converted',5,27)
% Will convert 27 time points starting at the 5th time point
% named converted_i0001.img, converted_i0002.img, ...,
% converted_i0027.img and their associated header and mat files.
%
% Assumes the data is stored as 003/I.001 through 003/I.999 and
% then 023/I.001 through 023/I.999, etc.
% Remember to skip the template images ie 1 time point for the
% first run.  You don't need to do this for subsequent runs.
%
% Modified to use spm write functions
% Modified to use spm_write_plane instead of spm_write_vol
%
% Souheil J. Inati  
% Dartmouth College
% September 2001
% souheil.inati@dartmouth.edu
%

if (nargin < 4)
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

% Is the first image the first or the last?
if (se_hdr.end_loc - se_hdr.start_loc) > 0
  scandir = 1;
else
  scandir = -1;
end

% Compute the M matrix
M = GE_createSPMmat(im_hdr,scandir);

% Create the template SPM volume structure
V.fname = '';
V.dim = [volSize spm_type('int16')]; % short ints
V.mat = M;
V.pinfo = [1 0 0]';

% Loop over the volumes until there is no more looping to be done.
for i = 1:nvols
  passnum = starttime + i - 1;

  % Read in the volume
  [imageVol, lastfile] = GE_readVolume(inDir, passnum, volSize, ...
				       pix_hdr.img_depth, im_offset);

  % Create analyze file
  vol_str = sprintf('000%d',i);
  vol_str = vol_str(length(vol_str)-3:length(vol_str));
  outName = [outStem sprintf('_i%s.img',vol_str)];
  V.fname = outName;
  V = spm_create_image(V);

  % Write out the SPM Volume
  % Don't use spm_write_vol because that sets a scale factor.
  for i = 1:Nz
    V = spm_write_plane(V,squeeze(imageVol(:,:,i)),i);
  end

  % Update to the screen
  fprintf('Wrote %s\n',outName);

end

% Done
fprintf('\nConversion Finished. \n');

return
