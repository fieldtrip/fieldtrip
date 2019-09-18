function hdr = load_analyze_hdr(hdrfile)
% hdr = load_analyze_hdr(hdrfile)
%
% hdrfile is the header file (eg, f001.hdr) or a basename (eg, f001).
%
% hdr.dime.dim(1) = number dims
% hdr.dime.dim(2) = number of items in first dim
% hdr.dime.dim(3) = number of items in second dim
% ...
%
% hdr.dime.datatype:
%  0    unk (???)
%  1    1 bit/pix
%  2    8 bits
%  4   16 bits
%  8   32 bits (signed int)
% 16   32 bits (floating pt)
% 32   64 bits (2 floats)
% 64   64 bits (double)
%
% hdr.dime.pixdim(1) = physical size of first dim (eg, 3.125 mm or 2000 ms)
% hdr.dime.pixdim(2) = ...
% 
% hdr.vox2ras is the vox2ras matrix if the hdrfile is accompanied
%   by a .mat file.
%


%
% load_analyze_hdr.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.6 $
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
  fprintf('hdr = load_analyze_hdr(hdrfile)\n');
  return;
end

% Try opening as big endian first
fp = fopen(hdrfile,'r','b');
if(fp == -1) 
  hdrfile0 = hdrfile;
  hdrfile = sprintf('%s.hdr',hdrfile);
  fp = fopen(hdrfile,'rb');
  if(fp == -1) 
    fprintf('ERROR: could not read %s or %s\n',hdrfile0,hdrfile);
    return;
  end
end


hdr.key.sizeof_hdr  = fread(fp,1,'int');
if(hdr.key.sizeof_hdr ~= 348)
  fclose(fp);
  % Now try opening as little endian
  fp = fopen(hdrfile,'r','l');
  hdr.key.sizeof_hdr  = fread(fp,1,'int');
  if(hdr.key.sizeof_hdr ~= 348)
    fclose(fp);
    fprintf('ERROR: %s: hdr size = %d, should be 348\n',...
	    hdrfile,hdr.key.sizeof_hdr);
    hdr = [];
    return;
  end
  hdr.endian = 'l';
else
  hdr.endian = 'b';
end

hdr.key.data_type     = fscanf(fp,'%c',10);

hdr.key.db_name       = fscanf(fp,'%c',18);
hdr.key.extents       = fread(fp,1,'int');
hdr.key.session_error = fread(fp,1,'short');
hdr.key.regular       = fscanf(fp,'%c',1);
hdr.key.hkey_un0      = fscanf(fp,'%c',1);

hdr.dime.dim        = fread(fp,8,'short');
hdr.dime.vox_units  = fscanf(fp,'%c',4);
hdr.dime.cal_units  = fscanf(fp,'%c',8);
hdr.dime.unused1    = fread(fp,1,'short');
hdr.dime.datatype   = fread(fp,1,'short');
hdr.dime.bitpix     = fread(fp,1,'short');
hdr.dime.dim_un0    = fread(fp,1,'short');
hdr.dime.pixdim     = fread(fp,8,'float');
hdr.dime.vox_offset = fread(fp,1,'float');
hdr.dime.roi_scale  = fread(fp,1,'float');
hdr.dime.funused1   = fread(fp,1,'float');
hdr.dime.funused2   = fread(fp,1,'float');
hdr.dime.cal_max    = fread(fp,1,'float');
hdr.dime.cal_min    = fread(fp,1,'float');
hdr.dime.compressed = fread(fp,1,'int');
hdr.dime.verified   = fread(fp,1,'int');
hdr.dime.glmax      = fread(fp,1,'int');
hdr.dime.glmin      = fread(fp,1,'int');

hdr.hist.descrip     = fscanf(fp,'%c',80);
hdr.hist.aux_file    = fscanf(fp,'%c',24);
hdr.hist.orient      = fscanf(fp,'%c',1);
hdr.hist.originator  = fscanf(fp,'%c',10); % Is this correct? Should be unsiged char
hdr.hist.generated   = fscanf(fp,'%c',10);
hdr.hist.scannum     = fscanf(fp,'%c',10);
hdr.hist.patient_id  = fscanf(fp,'%c',10);
hdr.hist.exp_date    = fscanf(fp,'%c',10);
hdr.hist.exp_time    = fscanf(fp,'%c',10);
hdr.hist.hist_un0    = fscanf(fp,'%c',3);
hdr.hist.views       = fread(fp,1,'int');
hdr.hist.vols_added  = fread(fp,1,'int');
hdr.hist.start_field = fread(fp,1,'int');
hdr.hist.field_skip  = fread(fp,1,'int');
hdr.hist.omax        = fread(fp,1,'int');
hdr.hist.omin        = fread(fp,1,'int');
hdr.hist.smax        = fread(fp,1,'int');
hdr.hist.smin        = fread(fp,1,'int');

fclose(fp);

% Load the .mat file if there
basename = hdrfile(1:end-4);
matfile = sprintf('%s.mat',basename);
fp = fopen(matfile,'r');
if(fp == -1) 
  hdr.vox2ras = diag([hdr.dime.pixdim(2:4)' 1]);
  % v = voxsize*ndim/2; -2 needed for 1-based
  v = abs(hdr.dime.pixdim .* (hdr.dime.dim+2)/2); 
  hdr.vox2ras(1,4) = +v(2); 
  hdr.vox2ras(2,4) = -v(3);
  hdr.vox2ras(3,4) = -v(4);
  hdr.have_vox2ras = 0;
else
  fclose(fp);
  tmp = load(matfile);
  hdr.vox2ras = tmp.M;
  hdr.have_vox2ras = 1;
end


return;





