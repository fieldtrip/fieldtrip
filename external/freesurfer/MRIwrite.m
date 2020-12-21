function err = MRIwrite(mri,fstring,datatype,permuteflag)
% err = MRIwrite(mri,fstring,datatype,permuteflag)
%
% Writes an MRI volume based on the fstring. fstring can be:
%  1. MGH file. Eg, f.mgh or f.mgz
%  2. BHDR file Eg, f.bhdr. Result will be written to a bfloat
%     volume, eg f_000.bfloat.
%  3. NIFIT file, Eg, f.nii, f.nii.gz (uncompressed and compressed)
%
% mri should be a structure like that read by MRIread.m The geometry
% (ie, direction cosines, voxel resolution, and P0 are all recomputed
% from mri.vox2ras0. So, if in the course of analysis, you changed
% mri.x_r, this change will not be reflected in the output volume.
%
% The only thing you need to fill-in in the mri struct is the mri.vol.
% All other fields will be filled in with defaulted values. 
% Fields are: vol, tr, vox2ras0, te, ti, flip_angle.
%
% 
% When writing in bhdr format, the default will be bfloat. If you want
% bshort, then set mri.outbext = 'bshort'. When a bhdr file is read in
% with MRIread(), mri.srcbext is set to either bshort or bfloat, so to
% keep the same precision set mri.outbext = mri.srcbext.  This only
% applies to bhdr format.
% 
% datatype can be uchar, short, int, float, double, ushort,
% uint. Only applies to nifti. Setting datatype to '' implies float.
%
% Note: you can write unpermuted data with permuteflag=0. If you leave
% out this flag or set it to anything non-zero, then it will be
% permuted. The permuteflag will not affect the vox2ras info. If you
% used permuteflag=0 when reading, then use the same when running
% MRIwrite().
% 

%
% MRIwrite.m
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


err = 1;

if(nargin < 2 | nargin > 4)
  fprintf('err = MRIwrite(mri,fstring,<datatype>,<permuteflag>)\n');
  return;
end
if(exist('datatype')~=1) datatype = ''; end  
if(isempty(datatype))   datatype = 'float'; end
if(exist('permuteflag')~=1) permuteflag = 1; end

if(~isfield(mri,'vol'))
  fprintf('ERROR: MRIwrite: structure does not have a vol field\n');
  return;
end
vsz = size(mri.vol);
nvsz = length(vsz);
if(nvsz ~= 4) vsz = [vsz ones(1,4-nvsz)]; end
if(~isfield(mri,'volsize'))
  mri.volsize = [vsz(1) vsz(2) vsz(3)];
end
if(~isfield(mri,'nframes'))  mri.nframes = vsz(4); end
if(~isfield(mri,'tr')) mri.tr = 0; end
if(~isfield(mri,'te')) mri.te = 0; end
if(~isfield(mri,'ti')) mri.ti = 0; end
if(~isfield(mri,'flip_angle')) mri.flip_angle = 0;end
if(~isfield(mri,'vox2ras0'))  mri.vox2ras0 = eye(4);end
if(~isfield(mri,'volres'))
    mri.volres = sqrt(sum(mri.vox2ras0(1:3,1:3).^2));
end
if(~isfield(mri,'scl_slope')) mri.scl_slope = 0; end
if(~isfield(mri,'scl_inter')) mri.scl_inter = 0; end
  
[fspec fstem fmt] = MRIfspec(fstring,0); % 0 = turn off checkdisk
if(isempty(fspec))
  fprintf('ERROR: could not determine format of %s\n',fstring);
  return;
end

switch(fmt)
 case {'mgh','mgz'} %----------- MGH/MGZ ------------%
  M = mri.vox2ras0;
  mr_parms = [mri.tr mri.flip_angle mri.te mri.ti];
  if(permuteflag)
    err = save_mgh(permute(mri.vol,[2 1 3 4]), fspec, M, mr_parms);  
  else
    err = save_mgh(mri.vol,fspec, M, mr_parms);  
  end
  return;
 case {'bhdr'} %----------- BHDR ------------%
  bmri.te = mri.te;
  bmri.tr = mri.tr;
  bmri.ti = mri.ti;
  bmri.flip_angle = mri.flip_angle;
  bmri.voldim = [size(mri.vol,1) size(mri.vol,2) size(mri.vol,3)];
  bmri.nframes = size(mri.vol,4);
  bmri.T = mri.vox2ras0;
  % Recompute voxel size based on vox2ras, to assure that
  % things only depend upon vox2ras0.
  xsize = sqrt(sum(mri.vox2ras0(:,1).^2)); 
  ysize = sqrt(sum(mri.vox2ras0(:,2).^2));
  zsize = sqrt(sum(mri.vox2ras0(:,3).^2));
  bmri.volres = [xsize ysize zsize];
  outbext = 'bfloat';
  if(isfield(mri,'outbext'))
    if(strcmp(mri.outbext,'bshort')) outbext = 'bshort'; end
  end
  err = fast_svbslice(mri.vol,fstem,[],outbext,bmri);
 
 case {'nii','nii.gz'} %----------- NIFTI! ------------%
  hdr.data_type       = '';
  hdr.db_name         = '';
  hdr.extents         = 0;
  hdr.session_error   = 0;
  hdr.regular         = '';
  hdr.dim_info        = '';
  
  % Note that the order is 2 1 3 4
  hdr.dim = [mri.volsize(2) mri.volsize(1) mri.volsize(3)  mri.nframes];
  hdr.intent_p1       = 0;
  hdr.intent_p2       = 0;
  hdr.intent_p3       = 0;
  hdr.intent_code     = 0;
  
  switch(datatype)
   case 'uchar',  hdr.datatype = 2;   hdr.bitpix = 8*1;
   case 'short',  hdr.datatype = 4;   hdr.bitpix = 8*2;
   case 'int',    hdr.datatype = 8;   hdr.bitpix = 8*4;
   case 'float',  hdr.datatype = 16;  hdr.bitpix = 8*4;
   case 'double', hdr.datatype = 64;  hdr.bitpix = 8*8;
   case 'char',   hdr.datatype = 256; hdr.bitpix = 8*1;
   case 'ushort', hdr.datatype = 512; hdr.bitpix = 8*2;
   case 'uint',   hdr.datatype = 768; hdr.bitpix = 8*4;
   otherwise,
    fprintf('ERROR: unrecognized data type %s\n',datatype);
    return;
  end
  hdr.slice_start     = 0;

  % volres is permuted in MRIread()
  hdr.pixdim          = [0 mri.volres([2 1 3]) mri.tr]; % physical units
  %hdr.pixdim          = [0 mri.volres mri.tr]; % physical units
  hdr.vox_offset      = 348; % will be set again
  hdr.scl_slope       = mri.scl_slope;
  hdr.scl_inter       = mri.scl_inter;
  
  hdr.slice_end       = 0;
  
  hdr.slice_code      = 0;
  hdr.xyzt_units = bitor(2,16); % 2=mm, 16=msec

  hdr.cal_max         = max(mri.vol(:));
  hdr.cal_min         = min(mri.vol(:));
  hdr.slice_duration  = 0;
  hdr.toffset         = 0;
  hdr.glmax           = 0;
  hdr.glmin           = 0;
  hdr.descrip         = sprintf('%-80s','FreeSurfer matlab');
  hdr.aux_file        = '';
  hdr.qform_code      = 1; % 1=NIFTI_XFORM_SCANNER_ANAT
  hdr.sform_code      = 1; % 1=NIFTI_XFORM_SCANNER_ANAT
  
  % Qform (must be 6dof)
  [b,c,d,x,y,z,qfac] = vox2rasToQform(mri.vox2ras0);
  hdr.pixdim(1)       = qfac;
  hdr.quatern_b       = b;
  hdr.quatern_c       = c;
  hdr.quatern_d       = d;
  hdr.quatern_x       = x;
  hdr.quatern_y       = y;
  hdr.quatern_z       = z;
  % Sform (can by any affine)
  hdr.srow_x          = mri.vox2ras0(1,:);
  hdr.srow_y          = mri.vox2ras0(2,:);
  hdr.srow_z          = mri.vox2ras0(3,:);
  hdr.intent_name     = 'huh?';
  hdr.magic           = 'n+1';

  % Note that the order is 2 1 3 4
  if(permuteflag)
    hdr.vol = permute(mri.vol,[2 1 3 4]);
  else
    hdr.vol = mri.vol;
  end
  err = save_nifti(hdr,fspec);
 otherwise
  fprintf('ERROR: format %s not supported\n',fmt);
  return;
end

if(err) fprintf('ERROR: saving %s \n',fstring); end

return;











