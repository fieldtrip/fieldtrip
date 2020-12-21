function mri = MRIread(fstring,headeronly,permuteflag)
% mri = MRIread(fstring,headeronly,permuteflag)
%
% Reads in a volume based on the fstring. fstring can be:
%  1. A stem, in which case the format and full file name is determined
%     by finding a file on disk called fstring.ext, where ext can be
%     either mgh, mgz, img, bhdr, nii, or nii.gz
%  2. MGH file. Eg, f.mgh or f.mgz
%  3. BVolume HDR file. Eg, f.bhdr 
%  4. Analyze, eg, f.img or f.hdr
%  5. NIFTI, eg, f.nii or f.nii.gz
%
% Creates a structure similar to the FreeSurfer MRI struct
% defined in mri.h. Times are in ms and angles are in radians.
% The vox2ras0 matrix is the matrix that converts a 0-based
% column, row, and slice to XYZ. vox2ras1 is the same with
% 1-based indices. The volume is rows, cols, slices frames,
% but the vox2ras expects col, row, slice.
%
% !!!! If you intend to use indices obtained from Matlab in FreeSurfer programs
%
% bear in mind that mri.vox(j+1, i+1, k+1) = mri(i,j,k)
%
% where mri(i,j,k) refers to indices as they are seen in scuba, tkmedit 
% or used in binaries such as mri_convert
%
% This happens because matlab uses row major (ie, the "fasted" dim
% goes from one row to the next), whereas C uses col major. 
% So if you simply load in a matrix and view it with imagesc, 
% it will appear to be transposed.
%
% Note: you can load unpermuted data with permuteflag=0. If you 
% leave out this flag or set it to anything non-zero, then 
% it will be permuted. The permuteflag will not affect the vox2ras
% info. If you used permuteflag=0, then use the same when 
% running MRIwrite().
%
% If headeronly=1, then the pixel data are not read in.
%
% If the input is a bhdr, then mri.srcbext is set to either bshort
% or bfloat, depending upon the precision of the input. If the
% input is not bhdr, then mri.srcbext will exist but be empty.
% See also MRIwrite() and mri.outbext.
%
% If the input is NIFTI, then mri.niftihdr is the nifti header
% If the input is ANALYZE, then mri.analyzehdr is the analyze header
%


%
% MRIread.m
%
% Original Author: Doug Greve
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


mri = [];

if(nargin < 1 | nargin > 3)
  fprintf('mri = MRIread(fstring,headeronly,permuteflag)\n');
  return;
end
if(exist('headeronly') ~=1) headeronly = 0; end
if(exist('permuteflag')~=1) permuteflag = 1; end

[fspec fstem fmt] = MRIfspec(fstring);
if(isempty(fspec))
  err = sprintf('ERROR: cannot determine format of %s (%s)\n',fstring,mfilename);
  error(err);
  return;
end

mri.srcbext = '';    % empty be default
mri.analyzehdr = []; % empty be default
mri.bhdr = []; % empty be default

%-------------- MGH ------------------------%
switch(fmt)
  case {'mgh','mgz'}
  [mri.vol, M, mr_parms, volsz] = load_mgh(fspec,[],[],headeronly);
  if(isempty(M))
    fprintf('ERROR: loading %s as MGH\n',fspec);
    mri = [];
    return;
  end
  if(~headeronly)
    if(permuteflag) mri.vol = permute(mri.vol,[2 1 3 4]); end
    volsz = size(mri.vol);
  else
    mri.vol = [];
    volsz(1:2) = [volsz(2) volsz(1)];
  end
  if(isempty(mr_parms)) mr_parms = zeros(4,1); end
  tr = mr_parms(1);
  flip_angle = mr_parms(2);
  te = mr_parms(3);
  ti = mr_parms(4);
%--------------- bshort/bfloat -----------------------%  
 case {'bhdr'}
  if(~headeronly)
    [mri.vol bmri] = fast_ldbslice(fstem);  
    if(isempty(mri.vol))
      fprintf('ERROR: loading %s as bvolume\n',fspec);
      mri = [];
      return;
    end
    volsz = size(mri.vol);
  else
    mri.vol = [];
    bmri = fast_ldbhdr(fstem);
    if(isempty(bmri))
      fprintf('ERROR: loading %s as bvolume\n',fspec);
      mri = [];
      return;
    end
    [nslices nrows ncols ntp] = fmri_bvoldim(fstem);
    volsz = [nrows ncols nslices ntp];
  end
  [nrows ncols ntp fs ns endian bext] = fmri_bfiledim(fstem);
  mri.srcbext = bext;
  M = bmri.T;
  tr = bmri.tr;
  flip_angle = bmri.flip_angle;
  te = bmri.te;
  ti = bmri.ti;
  mri.bhdr = bmri;
%------- analyze -------------------------------------   
 case {'img'}
  hdr = load_analyze(fspec,headeronly);
  if(isempty(hdr))
    fprintf('ERROR: loading %s as analyze\n',fspec);
    mri = [];
    return;
  end
  volsz = hdr.dime.dim(2:end);
  indnz = find(volsz~=0);
  volsz = volsz(indnz);
  volsz = volsz(:)'; % just make sure it's a row vect  
  if(~headeronly) 
    if(permuteflag) mri.vol = permute(hdr.vol,[2 1 3 4]); end
  else            
    mri.vol = [];
  end
  volsz([1 2]) = volsz([2 1]); % Make consistent. No effect when rows=cols
  tr = 1000*hdr.dime.pixdim(5); % msec
  flip_angle = 0;
  te = 0;
  ti = 0;
  hdr.vol = []; % already have it above, so clear it 
  M = vox2ras_1to0(hdr.vox2ras);
  mri.analyzehdr = hdr;
%------- nifti nii -------------------------------------   
 case {'nii','nii.gz'}
  hdr = load_nifti(fspec,headeronly);
  if(isempty(hdr))
    fprintf('ERROR: loading %s as analyze\n',fspec);
    mri = [];
    return;
  end
  volsz = hdr.dim(2:end);
  indnz = find(volsz~=0);
  volsz = volsz(indnz);
  volsz = volsz(:)'; % just make sure it's a row vect
  % This handles the case where data has > 4 dims
  % Just puts all data into dim 4.
  if(~headeronly) 
    hdr.vol = reshape(hdr.vol,[volsz(1) volsz(2) volsz(3) prod(volsz(4:end))]);
    if(permuteflag) 
      mri.vol = permute(hdr.vol,[2 1 3 4]); 
    else
      mri.vol = hdr.vol;
    end
  else mri.vol = [];
  end
  volsz([1 2]) = volsz([2 1]); % Make consistent. No effect when rows=cols
  tr = hdr.pixdim(5); % already msec
  flip_angle = 0;
  te = 0;
  ti = 0;
  hdr.vol = []; % already have it above, so clear it 
  M = hdr.vox2ras;
  mri.niftihdr = hdr;
%--------------------------------------------------- 
 otherwise
  fprintf('ERROR: format %s not supported\n',fmt);
  mri = [];
  return;
end
%--------------------------------------%

mri.fspec = fspec;
mri.pwd = pwd;

mri.flip_angle = flip_angle;
mri.tr  = tr;
mri.te  = te;
mri.ti  = ti;

% Assumes indices are 0-based. See vox2ras1 below for 1-based.  Note:
% MRIwrite() derives all geometry information (ie, direction cosines,
% voxel resolution, and P0 from vox2ras0. If you change other geometry
% elements of the structure, it will not be reflected in the output
% volume. Also note that vox2ras still maps Col-Row-Slice and not
% Row-Col-Slice.  Make sure to take this into account when indexing
% into matlab volumes (which are RCS). 
mri.vox2ras0 = M; 

% Dimensions not redundant when using header only
volsz(length(volsz)+1:4) = 1; % Make sure all dims are represented
mri.volsize = volsz(1:3); % only spatial components
mri.height  = volsz(1);   % Note that height (rows) is the first dimension
mri.width   = volsz(2);   % Note that width (cols) is the second dimension
mri.depth   = volsz(3);
mri.nframes = volsz(4);

%--------------------------------------------------------------------%
% Everything below is redundant in that they can be derivied from
% stuff above, but they are in the MRI struct defined in mri.h, so I
% thought I would add them here for completeness.  Note: MRIwrite()
% derives all geometry information (ie, direction cosines, voxel
% resolution, and P0 from vox2ras0. If you change other geometry
% elements below, it will not be reflected in the output volume.

mri.vox2ras = mri.vox2ras0;

mri.nvoxels = mri.height * mri.width * mri.depth; % number of spatial voxles
mri.xsize = sqrt(sum(M(:,1).^2)); % Col
mri.ysize = sqrt(sum(M(:,2).^2)); % Row
mri.zsize = sqrt(sum(M(:,3).^2)); % Slice

mri.x_r = M(1,1)/mri.xsize; % Col
mri.x_a = M(2,1)/mri.xsize;
mri.x_s = M(3,1)/mri.xsize;

mri.y_r = M(1,2)/mri.ysize; % Row
mri.y_a = M(2,2)/mri.ysize;
mri.y_s = M(3,2)/mri.ysize;

mri.z_r = M(1,3)/mri.zsize; % Slice
mri.z_a = M(2,3)/mri.zsize;
mri.z_s = M(3,3)/mri.zsize;

ic = [(mri.width)/2 (mri.height)/2 (mri.depth)/2 1]';
c = M*ic;
mri.c_r = c(1);
mri.c_a = c(2);
mri.c_s = c(3);
%--------------------------------------------------%

%-------- The stuff here is for convenience --------------

% 1-based vox2ras. Good for doing transforms in matlab
mri.vox2ras1 = vox2ras_0to1(M); 

% Matrix of direction cosines
mri.Mdc = [M(1:3,1)/mri.xsize M(1:3,2)/mri.ysize M(1:3,3)/mri.zsize];

% Vector of voxel resolutions (Row-Col-Slice)
mri.volres = [mri.ysize mri.xsize mri.zsize];

% Have to swap rows and columns back
voldim = [mri.volsize(2) mri.volsize(1) mri.volsize(3)]; %[ncols nrows nslices] 
volres = [mri.volres(2)  mri.volres(1)  mri.volres(3)];  %[dcol drow dslice] 
mri.tkrvox2ras = vox2ras_tkreg(voldim,volres);


return;











