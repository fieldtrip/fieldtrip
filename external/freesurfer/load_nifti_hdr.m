function hdr = load_nifti_hdr(niftifile)
% hdr = load_nifti_hdr(niftifile)
%
% Changes units to mm and msec.
% Creates hdr.sform and hdr.qform with the matrices in them.
% Creates hdr.vox2ras based on sform if valid, then qform.
% Does not and will not handle compressed. Compression is handled
% in load_nifti.m, which calls load_nifti_hdr.m after any
% decompression. 
%
% Endianness is returned as hdr.endian, which is either 'l' or 'b'. 
% When opening again, use fp = fopen(niftifile,'r',hdr.endian);
%
% Handles data structures with more than 32k cols by
% reading hdr.glmin = ncols when hdr.dim(2) < 0. This
% is FreeSurfer specific, for handling surfaces.
%


%
% load_nifti_hdr.m
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


hdr = [];

if(nargin ~= 1)
  fprintf('hdr = load_nifti_hdr(niftifile)\n');
  return;
end

% Try opening as big endian first
fp = fopen(niftifile,'r','b');
if(fp == -1) 
  fprintf('ERROR: could not read %s\n',niftifile);
  return;
end

hdr.sizeof_hdr  = fread(fp,1,'int');
if(hdr.sizeof_hdr ~= 348)
  fclose(fp);
  % Now try opening as little endian
  fp = fopen(niftifile,'r','l');
  hdr.sizeof_hdr  = fread(fp,1,'int');
  if(hdr.sizeof_hdr ~= 348)
    fclose(fp);
    fprintf('ERROR: %s: hdr size = %d, should be 348\n',...
	    niftifile,hdr.sizeof_hdr);
    hdr = [];
    return;
  end
  hdr.endian = 'l';
else
  hdr.endian = 'b';
end

hdr.data_type       = fscanf(fp,'%c',10);
hdr.db_name         = fscanf(fp,'%c',18);
hdr.extents         = fread(fp, 1,'int');
hdr.session_error   = fread(fp, 1,'short');
hdr.regular         = fread(fp, 1,'char');
hdr.dim_info        = fread(fp, 1,'char');
hdr.dim             = fread(fp, 8,'short');
hdr.intent_p1       = fread(fp, 1,'float');
hdr.intent_p2       = fread(fp, 1,'float');
hdr.intent_p3       = fread(fp, 1,'float');
hdr.intent_code     = fread(fp, 1,'short');
hdr.datatype        = fread(fp, 1,'short');
hdr.bitpix          = fread(fp, 1,'short');
hdr.slice_start     = fread(fp, 1,'short');
hdr.pixdim          = fread(fp, 8,'float'); % physical units
hdr.vox_offset      = fread(fp, 1,'float');
hdr.scl_slope       = fread(fp, 1,'float');
hdr.scl_inter       = fread(fp, 1,'float');
hdr.slice_end       = fread(fp, 1,'short');
hdr.slice_code      = fread(fp, 1,'char');
hdr.xyzt_units      = fread(fp, 1,'char');
hdr.cal_max         = fread(fp, 1,'float');
hdr.cal_min         = fread(fp, 1,'float');
hdr.slice_duration  = fread(fp, 1,'float');
hdr.toffset         = fread(fp, 1,'float');
hdr.glmax           = fread(fp, 1,'int');
hdr.glmin           = fread(fp, 1,'int');
hdr.descrip         = fscanf(fp,'%c',80);
hdr.aux_file        = fscanf(fp,'%c',24);
hdr.qform_code      = fread(fp, 1,'short');
hdr.sform_code      = fread(fp, 1,'short');
hdr.quatern_b       = fread(fp, 1,'float');
hdr.quatern_c       = fread(fp, 1,'float');
hdr.quatern_d       = fread(fp, 1,'float');
hdr.quatern_x       = fread(fp, 1,'float');
hdr.quatern_y       = fread(fp, 1,'float');
hdr.quatern_z       = fread(fp, 1,'float');
hdr.srow_x          = fread(fp, 4,'float');
hdr.srow_y          = fread(fp, 4,'float');
hdr.srow_z          = fread(fp, 4,'float');
hdr.intent_name     = fscanf(fp,'%c',16);
hdr.magic           = fscanf(fp,'%c',4);

fclose(fp);

% This is to accomodate structures with more than 32k cols
% FreeSurfer specific. See also mriio.c.
if(hdr.dim(2) < 0) 
  hdr.dim(2) = hdr.glmin; 
  hdr.glmin = 0; 
end

% look at xyz units and convert to mm if needed
xyzunits = bitand(hdr.xyzt_units,7); % 0x7
switch(xyzunits)
 case 1, xyzscale = 1000.000; % meters
 case 2, xyzscale =    1.000; % mm
 case 3, xyzscale =     .001; % microns
 otherwise, 
  fprintf('WARNING: xyz units code %d is unrecognized, assuming mm\n',xyzunits);
  xyzscale = 1; % just assume mm
end
hdr.pixdim(2:4) = hdr.pixdim(2:4) * xyzscale;
hdr.srow_x = hdr.srow_x * xyzscale;
hdr.srow_y = hdr.srow_y * xyzscale;
hdr.srow_z = hdr.srow_z * xyzscale;

% look at time units and convert to msec if needed
tunits = bitand(hdr.xyzt_units,3*16+8); % 0x38 
switch(tunits)
 case  8, tscale = 1000.000; % seconds
 case 16, tscale =    1.000; % msec
 case 32, tscale =     .001; % microsec
 otherwise,  tscale = 0; 
end
hdr.pixdim(5) = hdr.pixdim(5) * tscale;

% Change value in xyzt_units to reflect scale change
hdr.xyzt_units = bitor(2,16); % 2=mm, 16=msec

% Sform matrix
hdr.sform =  [hdr.srow_x'; 
	      hdr.srow_y'; 
	      hdr.srow_z';
	      0 0 0 1];

% Qform matrix - not quite sure how all this works,
% mainly just copied CH's code from mriio.c
b = hdr.quatern_b;
c = hdr.quatern_c;
d = hdr.quatern_d;
x = hdr.quatern_x;
y = hdr.quatern_y;
z = hdr.quatern_z;
a = 1.0 - (b*b + c*c + d*d);
if(abs(a) < 1.0e-7)
  a = 1.0 / sqrt(b*b + c*c + d*d);
  b = b*a;
  c = c*a;
  d = d*a;
  a = 0.0;
else
  a = sqrt(a);
end
r11 = a*a + b*b - c*c - d*d;
r12 = 2.0*b*c - 2.0*a*d;
r13 = 2.0*b*d + 2.0*a*c;
r21 = 2.0*b*c + 2.0*a*d;
r22 = a*a + c*c - b*b - d*d;
r23 = 2.0*c*d - 2.0*a*b;
r31 = 2.0*b*d - 2*a*c;
r32 = 2.0*c*d + 2*a*b;
r33 = a*a + d*d - c*c - b*b;
if(hdr.pixdim(1) < 0.0)
  r13 = -r13;
  r23 = -r23;
  r33 = -r33;
end
qMdc = [r11 r12 r13; r21 r22 r23; r31 r32 r33];
D = diag(hdr.pixdim(2:4));
P0 = [x y z]';
hdr.qform = [qMdc*D P0; 0 0 0 1];

if(hdr.sform_code ~= 0)
  % Use sform first
  hdr.vox2ras = hdr.sform;
elseif(hdr.qform_code ~= 0)
  % Then use qform first
  hdr.vox2ras = hdr.qform;
else
  fprintf('WARNING: neither sform or qform are valid in %s\n', ...
	  niftifile);
  D = diag(hdr.pixdim(2:4));
  P0 = [0 0 0]';
  hdr.vox2ras = [eye(3)*D P0; 0 0 0 1];
end

return;





